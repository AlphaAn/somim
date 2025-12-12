# Copyright 2025 Khac Duc An Thai
# This file is part of SOMIM.
# SOMIM is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# SOMIM is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with Foobar. If not, see <https://www.gnu.org/licenses/>.

from setuptools import setup, Extension
import pybind11
import sys

# Define the C++ extension module
ext_modules = [
    Extension(
        "somim_module",
        sorted([
            "somim_wrapper.cpp",      # The wrapper source
            "../core/libsomim.cpp",   # The core logic source
        ]),
        include_dirs=[
            pybind11.get_include(),   # Path to pybind11 headers
            "../core"                 # Path to libsomim.h
        ],
        language="c++",
        extra_compile_args=["-std=c++11", "-O3"],  # Optimization flags
    ),
]

setup(
    name="somim",
    version="2.0.0",
    description="Python interface for SOMIM (Search for Optimal Measurements by an Iterative Method))",
    ext_modules=ext_modules,
    zip_safe=False,
)