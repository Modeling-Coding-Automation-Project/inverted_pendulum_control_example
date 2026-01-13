#ifndef __FURUTA_PENDULUM_PID_CONTROLLER_HPP__
#define __FURUTA_PENDULUM_PID_CONTROLLER_HPP__

#include <cmath>

class FurutaPendulum_PID_Controller {
public:
  using FLOAT = double;

public:
  FurutaPendulum_PID_Controller(FLOAT Ts = 0.005, FLOAT theta_ref_rad = 0.0,
                                FLOAT alpha_ref_rad = 0.0, FLOAT v_limit = 12.0,
                                FLOAT alpha_ref_limit_rad = (12.0 * M_PI /
                                                             180.0),
                                FLOAT theta_to_alpha_sign = 1.0);

  ~FurutaPendulum_PID_Controller();

  void reset();

  void set_theta_reference_rad(FLOAT theta_ref_rad);
  void set_theta_reference_deg(FLOAT theta_ref_deg);

  void set_alpha_reference_rad(FLOAT alpha_ref_rad);
  void set_alpha_reference_deg(FLOAT alpha_ref_deg);

  // Compute motor voltage command [V]
  FLOAT calculate_manipulation(FLOAT theta, FLOAT alpha, FLOAT dtheta,
                               FLOAT dalpha);

public:
  // === Tuning gains ===
  // Theta PID output is interpreted as a commanded pendulum angle offset [rad].
  FLOAT kp_theta;
  FLOAT ki_theta;
  FLOAT kd_theta;

  // Pendulum PID outputs motor voltage [V] to track alpha_ref_cmd.
  FLOAT kp_alpha;
  FLOAT ki_alpha;
  FLOAT kd_alpha;

private:
  static FLOAT wrap_to_pi(FLOAT angle_rad);

  static bool should_integrate_u(FLOAT u_cmd, FLOAT u_cmd_sat, FLOAT err,
                                 FLOAT u_lim);

private:
  FLOAT _Ts;

  // References (can be changed later)
  FLOAT _theta_ref_rad;
  FLOAT _alpha_ref_rad;

  // Voltage saturation used for anti-windup
  FLOAT _v_limit;

  // Limit for commanded pendulum tilt generated from theta control.
  FLOAT _alpha_ref_limit_rad;

  // Mapping sign: if theta is positive, we want alpha_ref_cmd to be positive.
  FLOAT _theta_to_alpha_sign;

  // Integrator states
  FLOAT _int_theta;
  FLOAT _int_alpha;

  // Simple derivative filtering (on measurement)
  FLOAT _dalpha_filt;
  FLOAT _dalpha_tau;
  FLOAT _dtheta_filt;
  FLOAT _dtheta_tau;
};

#endif // __FURUTA_PENDULUM_PID_CONTROLLER_HPP__
