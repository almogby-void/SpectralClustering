from setuptools import Extension, setup

module = Extension('capi_project',
                  sources=[
                    #'spkmeans',
                    'functions.c',
                    'spkmeansmodule.c'
                  ])
setup(name='capi_project',
     version='1.0.0',
     description='Python wrapper for C project',
     ext_modules=[module])