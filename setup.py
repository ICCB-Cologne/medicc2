import sys
from pathlib import Path

import numpy as np
from Cython.Build import cythonize
from setuptools import Extension, setup

sys.path.append('fstlib/cext')

this_directory = Path(__file__).parent
print(this_directory)
print(this_directory.absolute())
long_description = (this_directory / "README.md").read_text()

setup(
    name='MEDICC2',
    version='0.5b1',
    author='Tom L Kaufmann, Roland F Schwarz, Marina Petkovic',
    author_email='tkau93@gmail.com, roland.f.schwarz@gmail.com, marina.55kovic@gmail.com',
    description='Minimum Event Distance for Intra-tumour Copy-number Comparisons',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://bitbucket.org/schwarzlab/medicc2',
    classifiers=[
                "Programming Language :: Python :: 3",
                "License :: OSI Approved :: MIT License",
                "Operating System :: OS Independent",
    ],
    packages=['medicc', 'fstlib', 'fstlib.cext'],
    scripts=['tools/medicc-convert-old-input.py',
             'tools/medicc-create-cn-fst.py',
             'medicc2.py',
             'tools/medicc-phase.py'],
    license='GPL-3',
    include_dirs=np.get_include(),
    install_requires=[
        'numpy>=1.20.1'
    ],
    package_data={
        "": ["objects/*"],
    },
    ext_modules = cythonize([
        Extension("fstlib.cext.pywrapfst", 
                  ["fstlib/cext/pywrapfst.pyx"],
                  include_dirs=['fstlib/cext'],
                  libraries=["fst", "fstfar", "fstscript", "fstfarscript"],
                  extra_compile_args=['-std=c++17'],
                  language = "c++"),

        Extension("fstlib.cext.ops", 
                  ["fstlib/cext/ops.pyx"],
                  include_dirs=['fstlib/cext'],
                  libraries=["fst", "fstfar", "fstscript", "fstfarscript"],
                  extra_compile_args=['-std=c++17'],
                  language = "c++")
    ])
)
