#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/iostream.h>
#include "QP.hpp"

namespace py = pybind11;

PYBIND11_MODULE(qp_module, m) {
    m.doc() = "Python bindings for custom QP solver";

    py::add_ostream_redirect(m, "cout_redirect");

    py::class_<QP>(m, "QP")
        .def(py::init<
              const Eigen::MatrixXd&,
              const Eigen::VectorXd&,
              const Eigen::MatrixXd&,
              const Eigen::VectorXd&,
              const Eigen::MatrixXd&,
              const Eigen::VectorXd&>(),
             py::arg("Q"),
             py::arg("q"),
             py::arg("G"),
             py::arg("h"),
             py::arg("A"),
             py::arg("b"))
        .def("solve", &QP::solve)
        .def("get_opt_value", &QP::get_opt_value)
        .def_readwrite("x", &QP::x)
        ;
}