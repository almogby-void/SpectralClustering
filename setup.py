from setuptools import Extension, setup

module = Extension('myproject',
                  sources=[
                    'spkmeans',
                    'spkmeansmodule.c'
                  ])
setup(name='myproject',
     version='1.0.0',
     description='Python wrapper for C project',
     ext_modules=[module])