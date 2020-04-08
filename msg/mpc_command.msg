# Fused local position in NED.

uint64 timestamp			# time since system start (microseconds)

# mpc solver flag
int32 flag

# mpc attitude
float32[4] q

# mpc control commands
float32 thrust
float32 tilt_angle
float32 torque_x
float32 torque_y
float32 torque_z
