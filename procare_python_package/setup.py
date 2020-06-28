# ----------------------------------------------------------------------------
# <                               ProCare                         	     >
# ----------------------------------------------------------------------------
# The MIT License (MIT)
#
# Copyright (c) 2020 Merveille Eguida
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
# IN THE SOFTWARE.
# ----------------------------------------------------------------------------


########################## OPEN3D ORIGINAL WORK ##############################
# ----------------------------------------------------------------------------
# -                        Open3D: www.open3d.org                            -
# ----------------------------------------------------------------------------
# The MIT License (MIT)
#
# Copyright (c) 2018 www.open3d.org
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
# IN THE SOFTWARE.
 

from setuptools import setup, find_packages

with open("README.md", "r") as f:
    long_description = f.read()

with open('requirements.txt', 'r') as f:
    lines = f.readlines()
install_requires = [line.strip() for line in lines if line != ('' or '\n')]


import procare.info as info
print(dir(info))
setup(
    name="procare", # Replace with your own username
    version=info.__version__,
    install_requires=install_requires,
    author="Merveille Eguida",
    author_email="keguida@unistra.fr",
    description="Point cloud registration for protein cavities",
    long_description=long_description,
    keywords="protein cavity point cloud pharmaophore comparison",
    long_description_content_type="text/markdown",
    url="https://github.com/kieguida/ProCare",
    packages=find_packages(),
    package_dir={'procare':'procare'},
    package_data={'procare':['open3d/*.so']},
    classifiers=[
	"Development Status :: 2 - Pre-Alpha",
	"Natural Language :: English",
	"Programming Language :: C",
        "Programming Language :: C++",
        "Programming Language :: Python :: 3",
	"Programming Language :: Python :: 3.6",
	"Programming Language :: Python :: 3.7",
	"Programming Language :: Python :: 3.8",
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX :: Linux",
	"Topic :: Scientific/Engineering",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        
    ],
    python_requires='>=3.6',
)
