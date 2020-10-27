""" Install PIPPy
"""
from distutils.core import setup

setup(name="PIPPy",
      version="0.1",
      packages=["pippy",
                "pippy.fitparser",
                "pippy.src",
                "examples"],
      package_dir={'pippy': "PIPPy"})#,
#      package_data={'pippy': ["tests/data/*.txt"]})
