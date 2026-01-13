#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "furuta_pendulum_pid_controller.hpp"

namespace furuta_pendulum_pid_controller_SIL {

namespace py = pybind11;

FurutaPendulum_PID_Controller controller;

void initialize(FurutaPendulum_PID_Controller::FLOAT Ts) {
  controller = FurutaPendulum_PID_Controller(Ts);
}

void set_theta_reference_rad(
    FurutaPendulum_PID_Controller::FLOAT theta_ref_rad) {

  controller.set_theta_reference_rad(theta_ref_rad);
}

FurutaPendulum_PID_Controller::FLOAT
calculate_manipulation(FurutaPendulum_PID_Controller::FLOAT theta,
                       FurutaPendulum_PID_Controller::FLOAT alpha,
                       FurutaPendulum_PID_Controller::FLOAT dtheta,
                       FurutaPendulum_PID_Controller::FLOAT dalpha) {

  return controller.calculate_manipulation(theta, alpha, dtheta, dalpha);
}

PYBIND11_MODULE(FurutaPendulumPidControllerSIL, m) {
  m.def("initialize", &initialize, "Initialize the module");
  m.def("set_theta_reference_rad", &set_theta_reference_rad,
        "Set theta reference [rad]");
  m.def("calculate_manipulation", &calculate_manipulation,
        "Compute motor voltage command [V]");
}

} // namespace furuta_pendulum_pid_controller_SIL
