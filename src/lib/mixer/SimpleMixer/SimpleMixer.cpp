/****************************************************************************
 *
 *   Copyright (c) 2012-2015 PX4 Development Team. All rights reserved.
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
 * @file mixer_simple.cpp
 *
 * Simple summing mixer.
 */

#include "SimpleMixer.hpp"

#include <stdio.h>
#include <stdlib.h>

#define debug(fmt, args...)	do { } while(0)
//#define debug(fmt, args...)	do { printf("[mixer] " fmt "\n", ##args); } while(0)

SimpleMixer::SimpleMixer(ControlCallback control_cb, uintptr_t cb_handle, mixer_simple_s *mixinfo) :
	Mixer(control_cb, cb_handle),
	_pinfo(mixinfo)
{
}

SimpleMixer::~SimpleMixer()
{
	if (_pinfo != nullptr) {
		free(_pinfo);
	}
}

unsigned SimpleMixer::set_trim(float trim)
{
	_pinfo->output_scaler.offset = trim;
	return 1;
}

unsigned SimpleMixer::get_trim(float *trim)
{
	*trim = _pinfo->output_scaler.offset;
	return 1;
}

int
SimpleMixer::parse_output_scaler(const char *buf, unsigned &buflen, mixer_scaler_s &scaler)
{
	int ret;
	int s[5];
	int n = -1;

	buf = findtag(buf, buflen, 'O');

	if ((buf == nullptr) || (buflen < 12)) {
		debug("output parser failed finding tag, ret: '%s'", buf);
		return -1;
	}

	if ((ret = sscanf(buf, "O: %d %d %d %d %d %n",
			  &s[0], &s[1], &s[2], &s[3], &s[4], &n)) != 5) {
		debug("out scaler parse failed on '%s' (got %d, consumed %d)", buf, ret, n);
		return -1;
	}

	buf = skipline(buf, buflen);

	if (buf == nullptr) {
		debug("no line ending, line is incomplete");
		return -1;
	}

	scaler.negative_scale	= s[0] / 10000.0f;
	scaler.positive_scale	= s[1] / 10000.0f;
	scaler.offset		= s[2] / 10000.0f;
	scaler.min_output	= s[3] / 10000.0f;
	scaler.max_output	= s[4] / 10000.0f;

	return 0;
}

int
SimpleMixer::parse_control_scaler(const char *buf, unsigned &buflen, mixer_scaler_s &scaler, uint8_t &control_group,
				  uint8_t &control_index)
{
	unsigned u[2];
	int s[5];

	buf = findtag(buf, buflen, 'S');

	if ((buf == nullptr) || (buflen < 16)) {
		debug("control parser failed finding tag, ret: '%s'", buf);
		return -1;
	}

	if (sscanf(buf, "S: %u %u %d %d %d %d %d",
		   &u[0], &u[1], &s[0], &s[1], &s[2], &s[3], &s[4]) != 7) {
		debug("control parse failed on '%s'", buf);
		return -1;
	}

	buf = skipline(buf, buflen);

	if (buf == nullptr) {
		debug("no line ending, line is incomplete");
		return -1;
	}

	control_group		= u[0];
	control_index		= u[1];
	scaler.negative_scale	= s[0] / 10000.0f;
	scaler.positive_scale	= s[1] / 10000.0f;
	scaler.offset		= s[2] / 10000.0f;
	scaler.min_output	= s[3] / 10000.0f;
	scaler.max_output	= s[4] / 10000.0f;

	return 0;
}

SimpleMixer *
SimpleMixer::from_text(Mixer::ControlCallback control_cb, uintptr_t cb_handle, const char *buf, unsigned &buflen)
{
	SimpleMixer *sm = nullptr;
	mixer_simple_s *mixinfo = nullptr;
	unsigned inputs;
	int used;
	const char *end = buf + buflen;
	char next_tag;

	/* enforce that the mixer ends with a new line */
	if (!string_well_formed(buf, buflen)) {
		return nullptr;
	}

	/* get the base info for the mixer */
	if (sscanf(buf, "M: %u%n", &inputs, &used) != 1) {
		debug("simple parse failed on '%s'", buf);
		goto out;
	}

	/* at least 1 input is required */
	if (inputs == 0) {
		debug("simple parse got 0 inputs");
		goto out;
	}

	buf = skipline(buf, buflen);

	if (buf == nullptr) {
		debug("no line ending, line is incomplete");
		goto out;
	}

	mixinfo = (mixer_simple_s *)malloc(MIXER_SIMPLE_SIZE(inputs));

	if (mixinfo == nullptr) {
		debug("could not allocate memory for mixer info");
		goto out;
	}

	mixinfo->control_count = inputs;

	/* find the next tag */
	next_tag = findnexttag(end - buflen, buflen);

	if (next_tag == 'S') {
		/* No output scalers specified. Use default values.
		 * Corresponds to:
		 * O:      10000  10000      0 -10000  10000
		 */
		mixinfo->output_scaler.negative_scale	= 1.0f;
		mixinfo->output_scaler.positive_scale	= 1.0f;
		mixinfo->output_scaler.offset		= 0.f;
		mixinfo->output_scaler.min_output	= -1.0f;
		mixinfo->output_scaler.max_output	= 1.0f;

	} else {
		if (parse_output_scaler(end - buflen, buflen, mixinfo->output_scaler)) {
			debug("simple mixer parser failed parsing out scaler tag, ret: '%s'", buf);
			goto out;
		}
	}

	for (unsigned i = 0; i < inputs; i++) {
		if (parse_control_scaler(end - buflen, buflen,
					 mixinfo->controls[i].scaler,
					 mixinfo->controls[i].control_group,
					 mixinfo->controls[i].control_index)) {
			debug("simple mixer parser failed parsing ctrl scaler tag, ret: '%s'", buf);
			goto out;
		}
	}

	sm = new SimpleMixer(control_cb, cb_handle, mixinfo);

	if (sm != nullptr) {
		mixinfo = nullptr;
		debug("loaded mixer with %d input(s)", inputs);

	} else {
		debug("could not allocate memory for mixer");
	}

out:

	if (mixinfo != nullptr) {
		free(mixinfo);
	}

	return sm;
}

unsigned
SimpleMixer::mix(float *outputs, unsigned space)
{

    // note: rest is assigned in multirotor mixer
    // load inputs
    float M,N, airspd;
    _control_cb(_cb_handle, 0, 1, M);
    M = math::constrain(M, -1.0f, 1.0f);

    _control_cb(_cb_handle, 0, 2, N);
    N = math::constrain(N, -1.0f, 1.0f);

    _control_cb(_cb_handle, 0, 5, airspd);
	airspd  = math::constrain( airspd, 1e-8f, 1.0f);

    // denormalize
    float M_MAX     = 2.0f;
    float AIRSPD_MAX= 40.0f;
    M       *= M_MAX;
    N       *= M_MAX;
    airspd  *= AIRSPD_MAX;

    float C_Me = 0.55604;
    float C_Nr = 0.055604;
    float S = 0.4266;
    float b = 2;
    float c_bar = 0.2;
    float delta_min = math::radians(-35.0f);
    float delta_max = math::radians(35.0f);

    float M_factor = C_Me * S * c_bar;
    float N_factor = C_Nr * S * b;

    float q_bar = 0.5f * 1.2f * airspd * airspd;
    float M_ = M_factor * q_bar;
    float N_ = N_factor * q_bar;


    // scale with airspeed to avoid bang-bang behaviour at low speeds
    float scale = math::constrain( (airspd - 4.0f)/6.0f, 0.0f, 1.0f);
    
    float delta_e = math::constrain( M/M_*scale, delta_min, delta_max );
    float delta_r = math::constrain( N/N_*scale,   delta_min, delta_max );

    outputs[0] = (2.0f * delta_e - (delta_max + delta_min))/(delta_max - delta_min);
    outputs[1] = (2.0f * delta_r - (delta_max + delta_min))/(delta_max - delta_min);
    outputs[2] = 0.0;

 return 3;
}

void
SimpleMixer::groups_required(uint32_t &groups)
{
	for (unsigned i = 0; i < _pinfo->control_count; i++) {
		groups |= 1 << _pinfo->controls[i].control_group;
	}
}

int
SimpleMixer::check()
{
	float junk;

	/* sanity that presumes that a mixer includes a control no more than once */
	/* max of 32 groups due to groups_required API */
	if (_pinfo->control_count > 32) {
		return -2;
	}

	/* validate the output scaler */
	int ret = scale_check(_pinfo->output_scaler);

	if (ret != 0) {
		return ret;
	}

	/* validate input scalers */
	for (unsigned i = 0; i < _pinfo->control_count; i++) {

		/* verify that we can fetch the control */
		if (_control_cb(_cb_handle,
				_pinfo->controls[i].control_group,
				_pinfo->controls[i].control_index,
				junk) != 0) {
			return -3;
		}

		/* validate the scaler */
		ret = scale_check(_pinfo->controls[i].scaler);

		if (ret != 0) {
			return (10 * i + ret);
		}
	}

	return 0;
}

float
SimpleMixer::scale(const mixer_scaler_s &scaler, float input)
{
	float output;

	if (input < 0.0f) {
		output = (input * scaler.negative_scale) + scaler.offset;

	} else {
		output = (input * scaler.positive_scale) + scaler.offset;
	}

	return math::constrain(output, scaler.min_output, scaler.max_output);
}

int
SimpleMixer::scale_check(mixer_scaler_s &scaler)
{
	if (scaler.offset > 1.001f) {
		return 1;
	}

	if (scaler.offset < -1.001f) {
		return 2;
	}

	if (scaler.min_output > scaler.max_output) {
		return 3;
	}

	if (scaler.min_output < -1.001f) {
		return 4;
	}

	if (scaler.max_output > 1.001f) {
		return 5;
	}

	return 0;
}
