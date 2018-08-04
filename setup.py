#!/usr/bin/env python  
import os  
import numpy
from distutils.core import setup, Extension  

MOD = '_count'  
setup(name=MOD, ext_modules=[Extension(MOD, ["./src/_count.cpp"], extra_compile_args=['-std=c++11'])], include_dirs=[numpy.get_include()])  
