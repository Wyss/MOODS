#!/bin/bash
# wrapper for python setup.py on OS X for x86_64 if build problems
# relate to architecture linking
# use:
#   x86_64_setup.sh build_ext --inplace
# is equivalent to: 
#   python setup.py build_ext --inplace
ARCHFLAGS='-arch x86_64' python setup.py $@