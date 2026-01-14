#include "furuta_pendulum_pid_controller.hpp"

FurutaPendulum_PID_Controller::FurutaPendulum_PID_Controller(
    FurutaPendulum_PID_Controller::FLOAT Ts,
    FurutaPendulum_PID_Controller::FLOAT theta_ref_rad,
    FurutaPendulum_PID_Controller::FLOAT alpha_ref_rad,
    FurutaPendulum_PID_Controller::FLOAT v_limit,
    FurutaPendulum_PID_Controller::FLOAT alpha_ref_limit_rad,
    FurutaPendulum_PID_Controller::FLOAT theta_to_alpha_sign)
    : kp_theta(FurutaPendulum_PID_Controller::KP_THETA_DEFAULT),
      ki_theta(FurutaPendulum_PID_Controller::KI_THETA_DEFAULT),
      kd_theta(FurutaPendulum_PID_Controller::KD_THETA_DEFAULT),
      kp_alpha(FurutaPendulum_PID_Controller::KP_ALPHA_DEFAULT),
      ki_alpha(FurutaPendulum_PID_Controller::KI_ALPHA_DEFAULT),
      kd_alpha(FurutaPendulum_PID_Controller::KD_ALPHA_DEFAULT),
      _Ts(static_cast<FLOAT>(Ts)),
      _theta_ref_rad(static_cast<FLOAT>(theta_ref_rad)),
      _alpha_ref_rad(static_cast<FLOAT>(alpha_ref_rad)),
      _v_limit(static_cast<FLOAT>(std::abs(v_limit))),
      _alpha_ref_limit_rad(static_cast<FLOAT>(std::abs(alpha_ref_limit_rad))),
      _theta_to_alpha_sign(theta_to_alpha_sign >= 0.0 ? 1.0 : -1.0),
      _int_theta(0.0), _int_alpha(0.0), _dalpha_filt(0.0),
      _dalpha_tau(FurutaPendulum_PID_Controller::DALPHA_TAU_DEFAULT),
      _dtheta_filt(0.0),
      _dtheta_tau(FurutaPendulum_PID_Controller::DTHETA_TAU_DEFAULT) {}

FurutaPendulum_PID_Controller::~FurutaPendulum_PID_Controller() {}

void FurutaPendulum_PID_Controller::reset() {
  this->_int_theta = 0.0;
  this->_int_alpha = 0.0;
  this->_dalpha_filt = 0.0;
  this->_dtheta_filt = 0.0;
}

void FurutaPendulum_PID_Controller::set_theta_reference_rad(
    FurutaPendulum_PID_Controller::FLOAT theta_ref_rad) {
  this->_theta_ref_rad = static_cast<FLOAT>(theta_ref_rad);
}

void FurutaPendulum_PID_Controller::set_theta_reference_deg(
    FurutaPendulum_PID_Controller::FLOAT theta_ref_deg) {
  this->_theta_ref_rad = static_cast<FLOAT>(theta_ref_deg) * M_PI / 180.0;
}

void FurutaPendulum_PID_Controller::set_alpha_reference_rad(
    FurutaPendulum_PID_Controller::FLOAT alpha_ref_rad) {
  this->_alpha_ref_rad = static_cast<FLOAT>(alpha_ref_rad);
}

void FurutaPendulum_PID_Controller::set_alpha_reference_deg(
    FurutaPendulum_PID_Controller::FLOAT alpha_ref_deg) {
  this->_alpha_ref_rad = static_cast<FLOAT>(alpha_ref_deg) * M_PI / 180.0;
}

FurutaPendulum_PID_Controller::FLOAT FurutaPendulum_PID_Controller::wrap_to_pi(
    FurutaPendulum_PID_Controller::FLOAT angle_rad) {
  const FLOAT two_pi = static_cast<FLOAT>(2.0 * M_PI);
  FLOAT x = static_cast<FLOAT>(std::fmod(angle_rad + M_PI, two_pi));
  if (x < 0.0) {
    x += two_pi;
  }
  return x - static_cast<FLOAT>(M_PI);
}

bool FurutaPendulum_PID_Controller::should_integrate_u(
    FurutaPendulum_PID_Controller::FLOAT u_cmd,
    FurutaPendulum_PID_Controller::FLOAT u_cmd_sat,
    FurutaPendulum_PID_Controller::FLOAT err,
    FurutaPendulum_PID_Controller::FLOAT u_lim) {
  if (u_cmd == u_cmd_sat) {
    return true;
  }
  const FLOAT eps = static_cast<FLOAT>(1e-12);
  if (u_cmd_sat >= u_lim - eps && err > 0.0) {
    return false;
  }
  if (u_cmd_sat <= -u_lim + eps && err < 0.0) {
    return false;
  }
  return true;
}

FurutaPendulum_PID_Controller::FLOAT
FurutaPendulum_PID_Controller::calculate_manipulation(
    FurutaPendulum_PID_Controller::FLOAT theta,
    FurutaPendulum_PID_Controller::FLOAT alpha,
    FurutaPendulum_PID_Controller::FLOAT dtheta,
    FurutaPendulum_PID_Controller::FLOAT dalpha) {
  FLOAT Ts = this->_Ts;
  if (Ts <= FurutaPendulum_PID_Controller::TS_MIN) {
    Ts = static_cast<FLOAT>(FurutaPendulum_PID_Controller::TS_MIN);
  }

  const FLOAT theta_meas = static_cast<FLOAT>(theta);
  const FLOAT alpha_meas = static_cast<FLOAT>(alpha);
  const FLOAT dtheta_meas = static_cast<FLOAT>(dtheta);
  const FLOAT dalpha_meas = static_cast<FLOAT>(dalpha);

  // --- Theta error (defined so theta>ref -> positive error) ---
  const FLOAT e_theta = wrap_to_pi(theta_meas - _theta_ref_rad);

  // --- Filtered derivatives (measurement) ---
  FLOAT dtheta_used = dtheta_meas;
  if (this->_dtheta_tau > 0.0) {
    const FLOAT a = this->_dtheta_tau / (this->_dtheta_tau + Ts);
    const FLOAT b = Ts / (this->_dtheta_tau + Ts);
    this->_dtheta_filt = a * this->_dtheta_filt + b * dtheta_meas;
    dtheta_used = this->_dtheta_filt;
  }

  FLOAT dalpha_used = dalpha_meas;
  if (this->_dalpha_tau > 0.0) {
    const FLOAT a = this->_dalpha_tau / (this->_dalpha_tau + Ts);
    const FLOAT b = Ts / (this->_dalpha_tau + Ts);
    this->_dalpha_filt = a * this->_dalpha_filt + b * dalpha_meas;
    dalpha_used = this->_dalpha_filt;
  }

  // --- Theta PID -> commanded pendulum angle (alpha_ref_cmd) ---
  const FLOAT alpha_offset_cmd =
      static_cast<FLOAT>(this->_theta_to_alpha_sign) *
      (kp_theta * e_theta + ki_theta * this->_int_theta +
       kd_theta * dtheta_used);

  const FLOAT alpha_ref_cmd =
      wrap_to_pi(this->_alpha_ref_rad + alpha_offset_cmd);

  const FLOAT alpha_ref_cmd_sat =
      std::max(-this->_alpha_ref_limit_rad,
               std::min(this->_alpha_ref_limit_rad, alpha_ref_cmd));

  // Pendulum error relative to commanded reference
  const FLOAT e_alpha = wrap_to_pi(alpha_ref_cmd_sat - alpha_meas);

  // --- Pendulum stabilization voltage contribution ---
  const FLOAT v_unsat =
      static_cast<FLOAT>((kp_alpha * e_alpha) + (ki_alpha * _int_alpha) -
                         (kd_alpha * dalpha_used));

  // --- Saturation ---
  const FLOAT v_sat = std::max(
      -this->_v_limit, std::min(this->_v_limit, static_cast<FLOAT>(v_unsat)));

  // --- Anti-windup (conditional integration) ---
  if (should_integrate_u(alpha_ref_cmd, alpha_ref_cmd_sat, e_theta,
                         this->_alpha_ref_limit_rad)) {
    this->_int_theta += e_theta * Ts;
  }

  if (should_integrate_u(v_unsat, v_sat, e_alpha, this->_v_limit)) {
    this->_int_alpha += e_alpha * Ts;
  }

  return v_sat;
}
