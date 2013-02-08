HYDRO-BIDE
==========

Microbial ecological constraint-based community models

The repository contains:
* 1 Hydrodynamic BIDE (birth,immigration,death,emigration) simulation model based on neutral theory.
It contains aspects of conventional chemostate models (e.g. Monod, Droop).
* 1 Biogeographical simulation-based Monod model.
* 1 Biogeogrpahical simulation-based Droop (Cell Quota) Model.
* example files (.png) of figures the script creates
* files containing metagenomic data.  
  cow rumen data obtained from MG-RAST  http://metagenomics.anl.gov/metagenomics.cgi?page=MetagenomeProject&project=24 
* Module holding plotting and analysis functions

License
=======

This code is available under a [BSD 2-Clause License](http://opensource.org/licenses/bsd-license.php).

Copyright (c) 2013, Kenneth Locey. All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Contact information
-------------------
Ken Locey's email: ken@weecology.org and locey@biology.usu.edu

*Go to [Ken's website](http://kenlocey.weecology.org)*


Programs and packages used by the authors
-------------------------------

Python version 2.6.5 or higher
Numpy
Matplotlib
Cython (optional)

Cython
------
Most scripts will eventually have associated Cython setup files.


Cython is a language that makes writing C extensions for the Python language as easy as Python itself.
It is based on the well-known Pyrex, but supports more cutting edge functionality and optimizations.
*http://www.cython.org/*  
setup tutorial: *http://docs.cython.org/src/userguide/tutorial.html*

Most any Python script can be run at C speeds by:
* replacing the .py extension in your script with .pyx
* create a setup file. To save the user a hassle, we have made several setup files for scripts in this repository.
* run this command from within the folder containing the .pyx file and the setup file:


        python setup_thescriptsname.py build_ext --inplace

iPython and iPython Notebook
----------------------------

We are currently rescripting copies of repository files to be used in iPython Notebook. Running scripts in
iPython notebook will enable users to read comments in a convenient form, run blocks of code, and visualize
results between blocks of code. 
