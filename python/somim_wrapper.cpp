/* Copyright 2025 Khac Duc An Thai
// This file is part of SOMIM.
SOMIM is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
SOMIM is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with Foobar. If not, see <https://www.gnu.org/licenses/>.
*/

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>         
#include <pybind11/complex.h>     
#include <pybind11/numpy.h>
#include "libsomim.h"             

namespace py = pybind11;
using namespace std;
typedef complex<double> Complex;

// Helper function to load the 3D NumPy array (E) into the C++ SOMIM object
void set_ensemble_from_numpy(SOMIM& somim, py::array_t<Complex> E_array) {
    // 1. Check array dimensions
    if (E_array.ndim() != 3) {
        throw std::runtime_error("Ensemble E must be a 3D array (N x N x J)");
    }
    auto E_info = E_array.request();
    // E_info.shape[0] is N, E_info.shape[1] is N, E_info.shape[2] is J
    
    // 2. Iterate and set Rho (E[n, m, j] -> somim->setRho(j, n, m, z))
    for (int j = 0; j < E_info.shape[2]; ++j) {
        for (int n = 0; n < E_info.shape[0]; ++n) {
            for (int m = 0; m < E_info.shape[1]; ++m) {
                // Access data pointer and calculate index
                size_t index = j * E_info.strides[2] / sizeof(Complex) + 
                               n * E_info.strides[0] / sizeof(Complex) + 
                               m * E_info.strides[1] / sizeof(Complex);

                Complex z = E_array.data()[index];
                somim.setRho(j, n, m, z);
            }
        }
    }
}


// PYBIND11 MODULE DEFINITION
PYBIND11_MODULE(somim_module, m) {
    m.doc() = "pybind11 wrapper for the SOMIM accessible information library";

    // 1. Bind the SOMIM class
    py::class_<SOMIM>(m, "SOMIM")
        // Constructor: SOMIM(n, j, k)
        .def(py::init<int, int, int>())

        // Bind core functions
        .def("getMI", &SOMIM::getMI, "Calculates the accessible information")
        .def("setDirectGrad", &SOMIM::setDirectGrad)
        .def("setTolerance", &SOMIM::setTolerance)

        // Add the utility function for setting the ensemble E from NumPy
        .def("set_ensemble_from_numpy", &set_ensemble_from_numpy,
             "Sets the N x N x J ensemble E from a NumPy array.")
             
        // Add a function to get the resulting POVM Q back as a NumPy array
        .def("get_povm_as_numpy", [](SOMIM &somim) {
            int N = somim.getN();
            int K = somim.getK(); // Final K after elimination
            
            // Create a 3D NumPy array (N x N x K)
            py::array_t<Complex> Q_array({N, N, K});
            auto Q_info = Q_array.request();
            Complex *Q_ptr = (Complex *)Q_info.ptr;

            for (int k = 0; k < K; ++k) {
                for (int n = 0; n < N; ++n) {
                    for (int m = 0; m < N; ++m) {
                        // Access data pointer and calculate index
                        size_t index = k * Q_info.strides[2] / sizeof(Complex) + 
                                       n * Q_info.strides[0] / sizeof(Complex) + 
                                       m * Q_info.strides[1] / sizeof(Complex);
                        
                        Q_ptr[index] = somim.getPOVM(k, n, m);
                    }
                }
            }
            return Q_array;
        }, "Retrieves the optimal POVM Q as an N x N x K NumPy array.");
}