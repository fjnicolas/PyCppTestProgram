#include <pybind11/pybind11.h>
#include "src/mygeek.h"  // Include your C++ class header

namespace py = pybind11;

PYBIND11_MODULE(pybind_wrapper, m) {
    py::class_<MyGeek>(m, "MyGeek")
        .def(py::init<int>())
        .def("HelloGeek", &MyGeek::HelloGeek);
}
