from setuptools import Extension, setup

module = Extension('spkmeansmodule',
                  sources=[
                    'spkmeans.c',
                    'spkmeansmodule.c'
                  ])
setup(name='spkmeansmodule',
     version='1.0.0',
     description='Python wrapper for C project',
     ext_modules=[module])