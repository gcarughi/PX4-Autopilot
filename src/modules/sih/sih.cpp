/****************************************************************************
 *
 *   Copyright (c) 2019 PX4 Development Team. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 * 3. Neither the name PX4 nor the names of its contributors may be
 *    used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
 * OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 * AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 * ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 ****************************************************************************/

/**
 * @file sih.cpp
 * Simulator in Hardware
 *
 * @author Romain Chiappinelli      <romain.chiap@gmail.com>
 *
 * Coriolis g Corporation - January 2019
 */

#include "sih.hpp"

#include <px4_platform_common/getopt.h>
#include <px4_platform_common/log.h>

#include <drivers/drv_pwm_output.h>         // to get PWM flags

#include <unistd.h>
#include <string.h>
#include <fcntl.h>
#include <termios.h>

using namespace math;
using namespace matrix;


int Sih::custom_command(int argc, char *argv[])
{
	return print_usage("unknown command");
}


int Sih::task_spawn(int argc, char *argv[])
{
	_task_id = px4_task_spawn_cmd("sih",
				      SCHED_DEFAULT,
				      SCHED_PRIORITY_MAX,
				      2600,
				      (px4_main_t)&run_trampoline,
				      (char *const *)argv);

	if (_task_id < 0) {
		_task_id = -1;
		return -errno;
	}

	return 0;
}

Sih *Sih::instantiate(int argc, char *argv[])
{

	Sih *instance = new Sih();

	if (instance == nullptr) {
		PX4_ERR("alloc failed");
	}

	return instance;
}

Sih::Sih() :
	ModuleParams(nullptr),
	_loop_perf(perf_alloc(PC_ELAPSED, "sih_execution")),
	_sampling_perf(perf_alloc(PC_ELAPSED, "sih_sampling"))
{
}

void Sih::run()
{
	// initialize parameters
	parameters_update_poll();

	init_variables();
	init_sensors();

	const hrt_abstime task_start = hrt_absolute_time();
	_last_run = task_start;
	_gps_time = task_start;
	_serial_time = task_start;

	px4_sem_init(&_data_semaphore, 0, 0);

	hrt_call_every(&_timer_call, LOOP_INTERVAL, LOOP_INTERVAL, timer_callback, &_data_semaphore);

	perf_begin(_sampling_perf);

	while (!should_exit()) {
		px4_sem_wait(&_data_semaphore);     // periodic real time wakeup

		perf_end(_sampling_perf);
		perf_begin(_sampling_perf);

		perf_begin(_loop_perf);

		inner_loop();   // main execution function

		perf_end(_loop_perf);
	}

	hrt_cancel(&_timer_call);   // close the periodic timer interruption
	px4_sem_destroy(&_data_semaphore);
}

// timer_callback() is used as a real time callback to post the semaphore
void Sih::timer_callback(void *sem)
{
	px4_sem_post((px4_sem_t *)sem);
}

// this is the main execution waken up periodically by the semaphore
void Sih::inner_loop()
{
	_now = hrt_absolute_time();
	_dt = (_now - _last_run) * 1e-6f;
	_last_run = _now;

	read_motors();

	//generate_force_and_torques();
	generate_force_and_torques_vtol();

	equations_of_motion();

	reconstruct_sensors_signals();

	send_IMU();

	if (_now - _gps_time >= 50000) { // gps published at 20Hz
		_gps_time = _now;
		send_gps();
	}

	// send uart message every 40 ms
	if (_now - _serial_time >= 40000) {
		_serial_time = _now;

		publish_sih();  // publish _sih message for debug purpose

		parameters_update_poll();   // update the parameters if needed
	}
}

void Sih::parameters_update_poll()
{
	// check for parameter updates
	if (_parameter_update_sub.updated()) {
		// clear update
		parameter_update_s pupdate;
		_parameter_update_sub.copy(&pupdate);

		// update parameters from storage
		updateParams();
		parameters_updated();
	}
}

// store the parameters in a more convenient form
void Sih::parameters_updated()
{
    _CHECK_ROBUST = (bool)_sih_robust.get();
    _DIST       = _sih_dist.get();

	_T_MAX      = _sih_t_max.get();
	_Q_MAX      = _sih_q_max.get();
	_L_ROLL     = _sih_l_roll.get();
	_L_PITCH    = _sih_l_pitch.get();
	_KDV        = _sih_kdv.get();
	_KDW        = _sih_kdw.get();
	_H0         = _sih_h0.get();

	_LAT0       = (double)_sih_lat0.get() * 1.0e-7;
	_LON0       = (double)_sih_lon0.get() * 1.0e-7;
	_COS_LAT0   = cosl(radians(_LAT0));

	_MASS       = _sih_mass.get();

	_G_I        = Vector3f(0.0f, 0.0f, _MASS * CONSTANTS_ONE_G);

	_I = diag(Vector3f(_sih_ixx.get(), _sih_iyy.get(), _sih_izz.get()));
	_I(0, 1) = _I(1, 0) = _sih_ixy.get();
	_I(0, 2) = _I(2, 0) = _sih_ixz.get();
	_I(1, 2) = _I(2, 1) = _sih_iyz.get();

	_I_inv = inv(_I);

	_mu_I = Vector3f(_sih_mu_x.get(), _sih_mu_y.get(), _sih_mu_z.get());
}

// initialization of the variables for the simulator
void Sih::init_variables()
{
	srand(1234);    // initialize the random seed once before calling generate_wgn()

	_p_I = Vector3f(0.0f, 0.0f, 0.0f);
	_v_I = Vector3f(0.0f, 0.0f, 0.0f);
	_v_B = Vector3f(0.0f, 0.0f, 0.0f);
	_q   = Quatf(1.0f, 0.0f, 0.0f, 0.0f);
	_w_B = Vector3f(0.0f, 0.0f, 0.0f);

    for( int i = 0; i<NB_MOTORS; i++ ){
        _u[i] = 0.0f;
    }

    for( int i = 0; i<8; i++ ){
        _u_main[i] = 0.0f;
    }

    for( int i = 0; i<6; i++ ){
        _u_aux[i] = 0.0f;
    }

    _Ft_B.setZero();
    _Mt_B.setZero();
    _Fa_B.setZero();
    _Ma_B.setZero();

    
    /*
     * vehicle data
     */
    // dimensions and rotor positions
    L_0    = 0.29f;
    l_1    = 0.1575f;
    l_3    = 0.105f;
    l_4    = 0.105f;
    h_0    = 0.015f;
    h_1    = 0.05f;
    b      = 2.0f;
    c_bar  = 0.2f;

    d_i(0,0) = -l_3;    d_i(0,1) = l_4;     d_i(0,2) = l_4;     d_i(0,3) = -l_3;
    d_i(1,0) = L_0;     d_i(1,1) = L_0;     d_i(1,2) = -L_0;     d_i(1,3) = -L_0;
    d_i(2,0) = -h_0;    d_i(2,1) = -h_0;    d_i(2,2) = -h_0;    d_i(2,3) = -h_0;

    d_ei(0,0) = -l_1;   d_ei(0,1) = l_1;    d_ei(0,2) = l_1;    d_ei(0,3) = -l_1;
    d_ei(1,0) = 0.0f;   d_ei(1,1) = 0.0f;   d_ei(1,2) = 0.0f;   d_ei(1,3) = 0.0f;
    d_ei(2,0) = -h_1;   d_ei(2,1) = -h_1;   d_ei(2,2) = -h_1;   d_ei(2,3) = -h_1;

    d_ri.setZero();

    rho = 1.2f;

    //actuator ranges
    chi_max     = math::radians(90.0f); 
    chi_min     = math::radians(-10.0f);
    delta_max   = math::radians(35.0f); 
    delta_min   = math::radians(-35.0f); 

    //aerodynamic surfaces
    S_W = 0.4266f;
    S_Fn= 0.011f;
    S_Fe= 0.055f;
    S_Fd= 0.055f;
    S_E = 0.0465f;
    S_R = 0.0744f;

    S_S = 0.4266f;

    //aerodynamic coefficients
    float coeff_W[8]{20.0f,0.227f,0.25f,5.62f,0.03f,0.2f,1.0f,0.025f};
    float coeff_F[8]{20.0f,0.0f,0.0f,0.0f,1.28f,0.0f,0.0f,1.28f};
    float coeff_E[8]{1.0f,0.698f,0.0f,0.885f,0.0f,1.2354f,0.3137f,0.0f};
    float coeff_R[8]{1.0f,0.698f,0.0f,0.885f,0.0f,1.2354f,0.3137f,0.0f};

    aero_coeff.col(0) = Vector<float,8>(coeff_W);
    aero_coeff.col(1) = Vector<float,8>(coeff_F);
    aero_coeff.col(2) = Vector<float,8>(coeff_E);
    aero_coeff.col(3) = Vector<float,8>(coeff_R);

    C_La = 0.058649f;
    C_Me = 0.55604f;
    C_Nr = 0.055604f;

    //position of aerodynamic surfaces
    aero_pos.col(0) = Vector3f( 0.0f, 0.5f, -0.015f );    // right wing
    aero_pos.col(1) = Vector3f( 0.0f, -0.5f, -0.015f );   // left wing
    aero_pos.col(2) = Vector3f( 0.036f, 0.0f, 0.025f );   // fuselage
    aero_pos.col(3) = Vector3f( -0.71f, 0.0f, -0.015f );  // elevator
    aero_pos.col(4) = Vector3f( -0.71f, 0.0f, -0.04f );   // rudder
}

void Sih::init_sensors()
{
	_vehicle_gps_pos.fix_type = 3;  // 3D fix
	_vehicle_gps_pos.satellites_used = 8;
	_vehicle_gps_pos.heading = NAN;
	_vehicle_gps_pos.heading_offset = NAN;
	_vehicle_gps_pos.s_variance_m_s = 0.5f;
	_vehicle_gps_pos.c_variance_rad = 0.1f;
	_vehicle_gps_pos.eph = 0.9f;
	_vehicle_gps_pos.epv = 1.78f;
	_vehicle_gps_pos.hdop = 0.7f;
	_vehicle_gps_pos.vdop = 1.1f;
}

// read the motor signals outputted from the mixer
void Sih::read_motors()
{
	actuator_outputs_s actuators_out;

	//if (_actuator_out_sub.update(&actuators_out)) {
	//	for (int i = 0; i < NB_MOTORS; i++) { // saturate the motor signals
	//		_u[i] = constrain((actuators_out.output[i] - PWM_DEFAULT_MIN) / (PWM_DEFAULT_MAX - PWM_DEFAULT_MIN), 0.0f, 1.0f);
	//	}
	//}

    // load main motor pwm
	if (_actuator_out_sub_main.update(&actuators_out)) {
		for (int i = 0; i < 8; i++) { // saturate the motor signals
			_u_main[i] = constrain((actuators_out.output[i] - PWM_DEFAULT_MIN) / (PWM_DEFAULT_MAX - PWM_DEFAULT_MIN), 0.0f, 1.0f);
		}

	}

    // load aux motor pwm
	if (_actuator_out_sub_aux.update(&actuators_out)) {
		for (int i = 0; i < 6; i++) { // saturate the motor signals
			_u_aux[i] = constrain((actuators_out.output[i] - PWM_DEFAULT_MIN) / (PWM_DEFAULT_MAX - PWM_DEFAULT_MIN), 0.0f, 1.0f);
		}
	}

}

// generate the motors thrust and torque in the body frame
void Sih::generate_force_and_torques()
{
    //_Ft_B = Thrust + ind. Drag 
    //_Mt_B = arm times (Thrust + ind. Drag ) + resisting torque + gyroscopic torque
    //matrix::Vector3f T_temp;
	_Ft_B = Vector3f(0.0f, 0.0f, -_T_MAX * (+_u[0] + _u[1] + _u[2] + _u[3]));
	_Mt_B = Vector3f(_L_ROLL * _T_MAX * (-_u[0] + _u[1] + _u[2] - _u[3]),
			 _L_PITCH * _T_MAX * (+_u[0] - _u[1] + _u[2] - _u[3]),
			 _Q_MAX * (+_u[0] + _u[1] - _u[2] - _u[3]));

	_Fa_B = -_KDV * _v_I;   // first order drag to slow down the aircraft
	_Ma_B = -_KDW * _w_B;   // first order angular damper

}

// generate the motors thrust and torque in the body frame
void Sih::generate_force_and_torques_vtol()
{
    /*
     * compute rotor forces and moments
     * note: _u[i] are PWM values in [0,1] 
     * we assume that the rotor thrusts are given as T[i] = _T_MAX * _u[i]
     * and the induced torques as Q[i] = _Q_MAX * _T_MAX * _u[i], where
     * _Q_MAX is specified as _Q_MAX = C_Q / C_T
     */

    Vector<float,4> T;
    //T(0) = _T_MAX * _u[0];
    //T(1) = _T_MAX * _u[1];
    //T(2) = _T_MAX * _u[2];
    //T(3) = _T_MAX * _u[3];
    T(0) = (powf(2.0f*_u_main[0] + 0.14676f, 2.0f) - 0.0821782f) / 0.355359f;
    T(1) = (powf(2.0f*_u_main[1] + 0.14676f, 2.0f) - 0.0821782f) / 0.355359f;
    T(2) = (powf(2.0f*_u_main[2] + 0.14676f, 2.0f) - 0.0821782f) / 0.355359f;
    T(3) = (powf(2.0f*_u_main[3] + 0.14676f, 2.0f) - 0.0821782f) / 0.355359f;

    float chi_l     = -( 2.0f *_u_main[4] - 1.7106f ) / 0.9602f;
    float chi_r     =  ( 2.0f *_u_main[5] - 0.2894f ) / 0.9602f;

    float delta_a   = _u_main[6] * (delta_min - delta_max) + delta_max;
    float delta_e   = _u_aux[0] * (delta_max - delta_min) + delta_min;
    float delta_r   = _u_aux[1] * (delta_max - delta_min) + delta_min;

    //printf("delta_a: %f\n",(double)delta_a);
    //printf("delta_e: %f\n",(double)delta_e);
    //printf("delta_r: %f\n",(double)delta_r);

    Dcmf R_chir(AxisAngle<float>(Vector3f(0,1,0), -chi_r));
    Dcmf R_chil(AxisAngle<float>(Vector3f(0,1,0), -chi_l));

    Dcmf R_chi;

    // forces and moments
    Vector3f _T;
    Vector3f _Q;
    _Ft_B.setZero();
    _Mt_B.setZero();

    for( int i=0; i<4; i++){
        if( i < 2){
            R_chi = R_chir;
        } else {
            R_chi = R_chil;
        }

        d_ri    = Vector3f( d_i.col(i) ) + R_chi * Vector3f( d_ei.col(i) );
        _T      = - T(i) * Vector3f( R_chi.col(2) );
        _Q      = _Q_MAX * T(i) * Vector3f( R_chi.col(2) );

        _Ft_B += _T;
        _Mt_B += (float)std::pow(-1,i+1)*_Q + d_ri.cross( _T );

    }
    
    /*
     * compute aerodynamic forces and moments
     */
    _Fa_B.setZero();
    _Ma_B.setZero();

    Vector3f v_a;
    Vector3f F_a;
    Vector3f M_a;
    float AoA;
    float c_d;
    float c_l;

    //right wing
    v_a = _v_B + _w_B.cross( aero_pos.col(0) );
    AoA = std::atan2( v_a(2), v_a(0) );
    compute_aero_coeff( &c_d, &c_l, 0, AoA );

    F_a = - 0.25f * rho * S_W * sqrtf( v_a(0) * v_a(0) + v_a(2) * v_a(2) ) *
            ( v_a * c_d + Vector3f( -v_a(2), 0.0f, v_a(0) ) * c_l );
    M_a = Vector3f( aero_pos.col(0) ).cross( F_a );

    _Fa_B += F_a;
    _Ma_B += M_a;

    //printf("F_a: %f %f %f\n",(double)F_a(0),(double)F_a(1),(double)F_a(2));
    //printf("M_a: %f %f %f\n",(double)M_a(0),(double)M_a(1),(double)M_a(2));

    //left wing
    v_a = _v_B + _w_B.cross( aero_pos.col(1) );
    AoA = std::atan2( v_a(2), v_a(0) );
    compute_aero_coeff( &c_d, &c_l, 0, AoA );

    F_a = - 0.25f * rho * S_W * sqrtf( v_a(0) * v_a(0) + v_a(2) * v_a(2) )
            * ( v_a * c_d + Vector3f( -v_a(2), 0.0f, v_a(0) ) * c_l );
    M_a = Vector3f( aero_pos.col(1) ).cross( F_a );

    _Fa_B += F_a;
    _Ma_B += M_a;
    //printf("F_a: %f %f %f\n",(double)F_a(0),(double)F_a(1),(double)F_a(2));
    //printf("M_a: %f %f %f\n",(double)M_a(0),(double)M_a(1),(double)M_a(2));

    //fuselage
    v_a = _v_B + _w_B.cross( aero_pos.col(2) );
    compute_aero_coeff( &c_d, &c_l, 1, 0.0f );

    F_a(0) = S_Fn * sqrtf( v_a(0) * v_a(0) ) * v_a(0) * c_d;
    F_a(1) = S_Fe * sqrtf( v_a(1) * v_a(1) ) * v_a(1) * c_d;
    F_a(2) = S_Fd * sqrtf( v_a(2) * v_a(2) ) * v_a(2) * c_d;

    F_a *= -0.5f * rho;
    M_a = Vector3f( aero_pos.col(2) ).cross( F_a );

    _Fa_B += F_a;
    _Ma_B += M_a;
    //printf("_v_I: %f %f %f\n",(double)_v_I(0),(double)_v_I(1),(double)_v_I(2));
    //printf("_v_B: %f %f %f\n",(double)_v_B(0),(double)_v_B(1),(double)_v_B(2));
    //printf("v_a: %f %f %f\n",(double)v_a(0),(double)v_a(1),(double)v_a(2));
    //printf("F_a: %f %f %f\n",(double)F_a(0),(double)F_a(1),(double)F_a(2));
    //printf("M_a: %f %f %f\n",(double)M_a(0),(double)M_a(1),(double)M_a(2));

    //elevator
    v_a = _v_B + _w_B.cross( aero_pos.col(3) );
    AoA = std::atan2( v_a(2), v_a(0) );
    compute_aero_coeff( &c_d, &c_l, 2, AoA );

    F_a = - 0.5f * rho * S_E * sqrtf( v_a(0) * v_a(0) + v_a(2) * v_a(2) )
            * ( v_a * c_d + Vector3f( -v_a(2), 0.0f, v_a(0) ) * c_l );
    M_a = Vector3f( aero_pos.col(3) ).cross( F_a );

    _Fa_B += F_a;
    _Ma_B += M_a;
    //printf("F_a: %f %f %f\n",(double)F_a(0),(double)F_a(1),(double)F_a(2));
    //printf("M_a: %f %f %f\n",(double)M_a(0),(double)M_a(1),(double)M_a(2));

    //rudder
    v_a = _v_B + _w_B.cross( aero_pos.col(4) );
    AoA = std::atan2( v_a(1), v_a(0) );
    compute_aero_coeff( &c_d, &c_l, 3, AoA );

    F_a = - 0.5f * rho * S_R * sqrtf( v_a(0) * v_a(0) + v_a(1) * v_a(1) )
            * ( v_a * c_d + Vector3f( -v_a(1), v_a(0), 0.0f ) * c_l );
    M_a = Vector3f( aero_pos.col(4) ).cross( F_a );


    _Fa_B += F_a;
    _Ma_B += M_a;
    //printf("F_a: %f %f %f\n",(double)F_a(0),(double)F_a(1),(double)F_a(2));
    //printf("M_a: %f %f %f\n",(double)M_a(0),(double)M_a(1),(double)M_a(2));
    //printf("\n");

    /*
     * control surfaces
     */
    float q_bar     = 0.5f * rho * _v_B.norm() * _v_B.norm();
    float L_delta   = C_La * S_S * b     * q_bar * delta_a;
    float M_delta   = C_Me * S_S * c_bar * q_bar * delta_e;
    float N_delta   = C_Nr * S_S * b     * q_bar * delta_r;

    _Ma_B += Vector3f( L_delta, M_delta, N_delta );

    //printf("M_delta: %f\n",(double)M_delta);
    
    _Ma_B(0) += _DIST;
    

    /*
     * modify values to test robustness
     */
    if( _CHECK_ROBUST ){
      for(int i=0; i<3; i++){
        float k_rand = (float)rand()/RAND_MAX + 0.5f;
        _Ft_B(i) = k_rand * _Ft_B(i);
        _Mt_B(i) = k_rand * _Mt_B(i);
        _Fa_B(i) = k_rand * _Fa_B(i);
        _Ma_B(i) = k_rand * _Ma_B(i);
      }
    }
//    printf("_Ft_B: %f %f %f\n",(double)_Ft_B(0),(double)_Ft_B(1),(double)_Ft_B(2));
//    printf("_Fa_B: %f %f %f\n",(double)_Fa_B(0),(double)_Fa_B(1),(double)_Fa_B(2));
//    printf("_Mt_B: %f %f %f\n",(double)_Mt_B(0),(double)_Mt_B(1),(double)_Mt_B(2));
//    printf("_Ma_B: %f %f %f\n",(double)_Ma_B(0),(double)_Ma_B(1),(double)_Ma_B(2));
//    printf("\n");

}

void Sih::compute_aero_coeff(float* C_D, float* C_L, const int surf, const float AoA)
{
    float k         = aero_coeff(0,surf);
    float a_stall   = aero_coeff(1,surf);
    float c_l0      = aero_coeff(2,surf);
    float c_la      = aero_coeff(3,surf);
    float c_d0      = aero_coeff(4,surf);
    float c_da      = aero_coeff(5,surf);
    float c_1       = aero_coeff(6,surf);
    float c_0       = aero_coeff(7,surf);

    // blending scale
    float num = 1.0f + (float)std::tanh( k * ( a_stall * a_stall - AoA * AoA ) ); 
    float den = 1.0f + (float)std::tanh( k * a_stall * a_stall );
    float sigma = num / den;

    // drag coefficient
    float c_d = sigma * ( c_d0 + c_da * AoA * AoA ) +     
                (1.0f - sigma) * ( c_0 + c_1 * sin( AoA ) * sin( AoA ) );

    // lift coefficient
    float c_l = sigma * ( c_l0 + c_la * AoA ) +     
                (1.0f - sigma) * ( c_1 * sin( 2.0f * AoA ) );

    *C_D = c_d;
    *C_L = c_l;
}

// apply the equations of motion of a rigid body and integrate one step
void Sih::equations_of_motion()
{
	_C_IB = _q.to_dcm(); // body to inertial transformation

	// Equations of motion of a rigid body
	_p_I_dot = _v_I;                        // position differential
	_v_I_dot = (_G_I + _C_IB * (_Fa_B + _Ft_B)) / _MASS;   // conservation of linear momentum
	_q_dot = _q.derivative1(_w_B);              // attitude differential
	_w_B_dot = _I_inv * (_Mt_B + _Ma_B - _w_B.cross(_I * _w_B)); // conservation of angular momentum
    //printf("_v_I_dot: %f %f %f\n",(double)_v_I_dot(0),(double)_v_I_dot(1),(double)_v_I_dot(2));


	// fake ground, avoid free fall
	if (_p_I(2) > 0.0f && (_v_I_dot(2) > 0.0f || _v_I(2) > 0.0f)) {
		if (!_grounded) {    // if we just hit the floor
			// for the accelerometer, compute the acceleration that will stop the vehicle in one time step
			_v_I_dot = -_v_I / _dt;

		} else {
			_v_I_dot.setZero();
		}

		_v_I.setZero();
		_w_B.setZero();

		_grounded = true;

	} else {
		// integration: Euler forward
		_p_I = _p_I + _p_I_dot * _dt;
		_v_I = _v_I + _v_I_dot * _dt;
		_q = _q + _q_dot * _dt; // as given in attitude_estimator_q_main.cpp
		_q.normalize();
		_w_B = _w_B + _w_B_dot * _dt;
		_grounded = false;
	}

    _v_B = _C_IB.transpose() * _v_I;

    //printf("_p_I: %f %f %f\n",(double)_p_I(0),(double)_p_I(1),(double)_p_I(2));
    //printf("_v_I: %f %f %f\n",(double)_v_I(0),(double)_v_I(1),(double)_v_I(2));
    //printf("_q: %f %f %f %f\n",(double)_q(0),(double)_q(1),(double)_q(2),(double)_q(3));
    //printf("_w_B: %f %f %f\n",(double)_w_B(0),(double)_w_B(1),(double)_w_B(2));
    //printf("\n");
}

// reconstruct the noisy sensor signals
void Sih::reconstruct_sensors_signals()
{
	// The sensor signals reconstruction and noise levels are from
	// Bulka, Eitan, and Meyer Nahon. "Autonomous fixed-wing aerobatics: from theory to flight."
	// In 2018 IEEE International Conference on Robotics and Automation (ICRA), pp. 6573-6580. IEEE, 2018.

	// IMU
	_acc = _C_IB.transpose() * (_v_I_dot - Vector3f(0.0f, 0.0f, CONSTANTS_ONE_G)) + noiseGauss3f(0.5f, 1.7f, 1.4f);
	_gyro = _w_B + noiseGauss3f(0.14f, 0.07f, 0.03f);
	_mag = _C_IB.transpose() * _mu_I + noiseGauss3f(0.02f, 0.02f, 0.03f);

	// barometer
	float altitude = (_H0 - _p_I(2)) + generate_wgn() * 0.14f; // altitude with noise
	_baro_p_mBar = CONSTANTS_STD_PRESSURE_MBAR *        // reconstructed pressure in mBar
		       powf((1.0f + altitude * TEMP_GRADIENT / T1_K), -CONSTANTS_ONE_G / (TEMP_GRADIENT * CONSTANTS_AIR_GAS_CONST));
	_baro_temp_c = T1_K + CONSTANTS_ABSOLUTE_NULL_CELSIUS + TEMP_GRADIENT * altitude; // reconstructed temperture in celcius

	// GPS
	_gps_lat_noiseless = _LAT0 + degrees((double)_p_I(0) / CONSTANTS_RADIUS_OF_EARTH);
	_gps_lon_noiseless = _LON0 + degrees((double)_p_I(1) / CONSTANTS_RADIUS_OF_EARTH) / _COS_LAT0;
	_gps_alt_noiseless = _H0 - _p_I(2);

	_gps_lat = _gps_lat_noiseless + (double)(generate_wgn() * 7.2e-6f); // latitude in degrees
	_gps_lon = _gps_lon_noiseless + (double)(generate_wgn() * 1.75e-5f); // longitude in degrees
	_gps_alt = _gps_alt_noiseless + generate_wgn() * 1.78f;
	_gps_vel = _v_I + noiseGauss3f(0.06f, 0.077f, 0.158f);
}

void Sih::send_IMU()
{
	// gyro
	{
		static constexpr float scaling = 1000.0f;
		_px4_gyro.set_scale(1 / scaling);
		_px4_gyro.set_temperature(T1_C);
		_px4_gyro.update(_now, _gyro(0) * scaling, _gyro(1) * scaling, _gyro(2) * scaling);
	}

	// accel
	{
		static constexpr float scaling = 1000.0f;
		_px4_accel.set_scale(1 / scaling);
		_px4_accel.set_temperature(T1_C);
		_px4_accel.update(_now, _acc(0) * scaling, _acc(1) * scaling, _acc(2) * scaling);
	}

	// magnetometer
	{
		static constexpr float scaling = 1000.0f;
		_px4_mag.set_scale(1 / scaling);
		_px4_mag.set_temperature(T1_C);
		_px4_mag.update(_now, _mag(0) * scaling, _mag(1) * scaling, _mag(2) * scaling);
	}

	// baro
	{
		_px4_baro.set_temperature(_baro_temp_c);
		_px4_baro.update(_now, _baro_p_mBar);
	}
}

void Sih::send_gps()
{
	_vehicle_gps_pos.timestamp = _now;
	_vehicle_gps_pos.lat = (int32_t)(_gps_lat * 1e7);       // Latitude in 1E-7 degrees
	_vehicle_gps_pos.lon = (int32_t)(_gps_lon * 1e7); // Longitude in 1E-7 degrees
	_vehicle_gps_pos.alt = (int32_t)(_gps_alt * 1000.0f); // Altitude in 1E-3 meters above MSL, (millimetres)
	_vehicle_gps_pos.alt_ellipsoid = (int32_t)(_gps_alt * 1000); // Altitude in 1E-3 meters bove Ellipsoid, (millimetres)
	_vehicle_gps_pos.vel_ned_valid = true;              // True if NED velocity is valid
	_vehicle_gps_pos.vel_m_s = sqrtf(_gps_vel(0) * _gps_vel(0) + _gps_vel(1) * _gps_vel(
			1)); // GPS ground speed, (metres/sec)
	_vehicle_gps_pos.vel_n_m_s = _gps_vel(0);           // GPS North velocity, (metres/sec)
	_vehicle_gps_pos.vel_e_m_s = _gps_vel(1);           // GPS East velocity, (metres/sec)
	_vehicle_gps_pos.vel_d_m_s = _gps_vel(2);           // GPS Down velocity, (metres/sec)
	_vehicle_gps_pos.cog_rad = atan2(_gps_vel(1),
					 _gps_vel(0)); // Course over ground (NOT heading, but direction of movement), -PI..PI, (radians)

	_vehicle_gps_pos_pub.publish(_vehicle_gps_pos);
}

void Sih::publish_sih()
{
	// publish angular velocity groundtruth
	_vehicle_angular_velocity_gt.timestamp = hrt_absolute_time();
	_vehicle_angular_velocity_gt.xyz[0] = _w_B(0); // rollspeed;
	_vehicle_angular_velocity_gt.xyz[1] = _w_B(1); // pitchspeed;
	_vehicle_angular_velocity_gt.xyz[2] = _w_B(2); // yawspeed;

	_vehicle_angular_velocity_gt_pub.publish(_vehicle_angular_velocity_gt);

	// publish attitude groundtruth
	_att_gt.timestamp = hrt_absolute_time();
	_att_gt.q[0] = _q(0);
	_att_gt.q[1] = _q(1);
	_att_gt.q[2] = _q(2);
	_att_gt.q[3] = _q(3);

	_att_gt_pub.publish(_att_gt);

	_gpos_gt.timestamp = hrt_absolute_time();
	_gpos_gt.lat = _gps_lat_noiseless;
	_gpos_gt.lon = _gps_lon_noiseless;
	_gpos_gt.alt = _gps_alt_noiseless;
	_gpos_gt.vel_n = _v_I(0);
	_gpos_gt.vel_e = _v_I(1);
	_gpos_gt.vel_d = _v_I(2);

	_gpos_gt_pub.publish(_gpos_gt);
}

float Sih::generate_wgn()   // generate white Gaussian noise sample with std=1
{
	// algorithm 1:
	// float temp=((float)(rand()+1))/(((float)RAND_MAX+1.0f));
	// return sqrtf(-2.0f*logf(temp))*cosf(2.0f*M_PI_F*rand()/RAND_MAX);
	// algorithm 2: from BlockRandGauss.hpp
	static float V1, V2, S;
	static bool phase = true;
	float X;

	if (phase) {
		do {
			float U1 = (float)rand() / RAND_MAX;
			float U2 = (float)rand() / RAND_MAX;
			V1 = 2.0f * U1 - 1.0f;
			V2 = 2.0f * U2 - 1.0f;
			S = V1 * V1 + V2 * V2;
		} while (S >= 1.0f || fabsf(S) < 1e-8f);

		X = V1 * float(sqrtf(-2.0f * float(logf(S)) / S));

	} else {
		X = V2 * float(sqrtf(-2.0f * float(logf(S)) / S));
	}

	phase = !phase;
	return X;
}

// generate white Gaussian noise sample vector with specified std
Vector3f Sih::noiseGauss3f(float stdx, float stdy, float stdz)
{
	return Vector3f(generate_wgn() * stdx, generate_wgn() * stdy, generate_wgn() * stdz);
}

int sih_main(int argc, char *argv[])
{
	return Sih::main(argc, argv);
}

int Sih::print_usage(const char *reason)
{
	if (reason) {
		PX4_WARN("%s\n", reason);
	}

	PRINT_MODULE_DESCRIPTION(
		R"DESCR_STR(
### Description
This module provide a simulator for quadrotors running fully
inside the hardware autopilot.

This simulator subscribes to "actuator_outputs" which are the actuator pwm
signals given by the mixer.

This simulator publishes the sensors signals corrupted with realistic noise
in order to incorporate the state estimator in the loop.

### Implementation
The simulator implements the equations of motion using matrix algebra.
Quaternion representation is used for the attitude.
Forward Euler is used for integration.
Most of the variables are declared global in the .hpp file to avoid stack overflow.


)DESCR_STR");

    PRINT_MODULE_USAGE_NAME("sih", "simulation");
    PRINT_MODULE_USAGE_COMMAND("start");
    PRINT_MODULE_USAGE_DEFAULT_COMMANDS();

    return 0;
}
