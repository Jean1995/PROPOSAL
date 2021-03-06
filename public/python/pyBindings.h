#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#define DEF_PY_PRINT(cls) .def("__str__", &py_print<cls>)

namespace py = pybind11;

template <class T>
std::string py_print(const T& t) {
    std::stringstream ss;
    ss << t;
    return ss.str();
}
