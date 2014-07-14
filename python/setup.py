"""
When installing on OS X, you must build python with the same compiler 
as you build libpssm.a with, so basically, use the system python or 
build python from scratch from source or with from source with Homebrew
otherwise you may get unknown symbol errors
"""
from distutils.core import setup, Extension
from os.path import abspath, dirname, join

python_path = abspath(dirname(__file__))
moods_path=abspath(join(dirname(python_path), 'src'))

module1 = Extension('MOODS._cmodule',
                    sources = ['cport.cc'],
                    include_dirs=[moods_path],
                    library_dirs=[moods_path],
                    libraries=['pssm']
                    )

module2 = Extension('MOODS.moody',
                    sources = ['moods-py.cc', join(moods_path, 'mlf.cpp')],
                    include_dirs=[moods_path],
                    library_dirs = [moods_path],
                    libraries=['pssm']
                    )

import sys
import platform
if sys.platform == 'darwin' and 'clang' in platform.python_compiler().lower():
    print("Patching for Clang")  
    from distutils.sysconfig import get_config_vars
    import re
    res = get_config_vars()
    for key, value in res.items():
        if isinstance(value, str):
            if '-bundle' in value:
                flags = value
                flags = re.sub('-bundle', '-dynamiclib', flags)
                res[key] = flags
                print("XXX", key, ':', value, '-->', flags)
    print("End patching for Clang") 

setup (name = 'MOODS',
       version = '1.0',
       description = 'This is a demo package',
       packages = ['MOODS',],
       ext_modules = [module1, module2]
       )