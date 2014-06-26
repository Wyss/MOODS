from distutils.core import setup, Extension


module1 = Extension('MOODS._cmodule',
                    sources = ['cport.cc'],
                    include_dirs=["../src"],
                    library_dirs=["../src"],
                    libraries=['pssm'])

module2 = Extension('MOODS.moody',
                    sources = ['moods-py.c', '../src/pssm_algorithms.c'],
                    include_dirs=['../src', '../klib']
                    )

setup (name = 'MOODS',
       version = '1.0',
       description = 'This is a demo package',
       packages = ['MOODS',],
       ext_modules = [module1, module2])