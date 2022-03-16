from setuptools import setup, Extension

setup(
    name='myproject',
    version='1.0.0',
    description="Project",
    ext_modules=[Extension('myproject', sources=['spkmeansmodule.c'])])
