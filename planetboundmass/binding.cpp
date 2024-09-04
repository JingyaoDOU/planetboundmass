#include "bound.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(boundcpp, m) {
  py::class_<Boundcpp>(m, "Boundcpp")
      .def(py::init<>())
      .def("boundcpp", &Boundcpp::find_bound);
}
