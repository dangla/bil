Program
=======

Bil is a modeling platform based on finite volume/element methods. Bil is distributed under the terms of the GNU General Public License (GnuGPL). Bil can be downloaded from the URL: <https://github.com/ifsttar/bil>.


Building requirements
=====================

Building Bil from source requires:

  - Make (<http://www.gnu.org/software/make>)  
  - C and C++ compilers (<http://gcc.gnu.org>)  
  - Fortran compiler (<http://gcc.gnu.org/fortran>)
  
Building Bil on system other than Linux-based OS requires in addition

  - CMake (<http://www.cmake.org>)

Building the documentation from source requires:

  - Texinfo (<http://www.gnu.org/software/texinfo>)
  
and optionally

  - Doxygen (<http://www.stack.nl/~dimitri/doxygen>)
  

List of folders
===============

        ./base      reference solution and data bases
        ./bin       binaries  
        ./cmake     macros and functions for cmake  
        ./doc       documentations  
        ./examples  examples of input data files  
        ./src       sources  
        ./lib       libraries  
        ./scripts   some utility scripts  


Building and Installation
=========================

There are two methods to build Bil from the command line:

  - from make (native method suited for linux-based OS)
  - from cmake (suited for any OS like Windows or MacOS)
  

1. Build Bil from make (native method)
--------------------------------------
To build Bil and/or the documentation, use the following commands from the Bil's source directory:

        make      (build the binary file bil and the documentation)  
        make bin  (build the binary file only)  
        make doc  (build the documentation only)  

The binary file is created in the folder `./bin`, the documentation (files info, ps, pdf, txt, html) is created in `./doc`.

To install Bil, use the following commands (requires root permission):

        sudo make install      (install the binary and the documentation)  
        sudo make install-bin  (install the binary only)  
        sudo make install-doc  (install the documentation only)  

By default the files are installed in PREFIX = /usr/local

        bil       PREFIX/bin  
        bil.1     PREFIX/man/man1  
        bil.info  PREFIX/info  

PREFIX is defined in "Makefile". Unless essential it is not recommended to change the location PREFIX. If however you want another location, change the variable PREFIX in "Makefile" and  run "make" again. If you don't have root permissions use an alias instead (i.e. enter anywhere: alias bil='absolutepathtobilfolder/bin/bil-I.J').

Once the installation is completed, running 

        make clean  

will delete all the local files previously created.
  
  
2. Build Bil from cmake (command line)
--------------------------------------

Create a build directory, for example as a subdirectory of Bil's source directory:

        mkdir build  

Run cmake from within the build directory, pointing to Bil's source directory:

        cd build  
        cmake ..  

To build Bil then simply type:

        make  

the executable is copied in ../bin and the shared library in ../lib.

To install Bil type (this may require root permission):

        make install  

To change build options you can specify options directly on the command line, for example

        cd build  
        cmake -DCMAKE_INSTALL_PREFIX=/opt ..  

will change the installation directory (by default it is /usr/local on linux).

You can configure a debug build with

        cd build  
        cmake -DCMAKE_BUILD_TYPE=Debug ..  
        make  
        make install  
    
Note that the CMAKE_BUILD_TYPE variable is saved in the cache (CMakeCache.txt). On subsequent calls to cmake CMAKE_BUILD_TYPE is set to its cache value or to "Release"  the first time cmake is called.

You can keep multiple builds with different build options at the same time. Below the directory "build" you can place as many target directories for out-of-source build modes as you want. For example
  
        build/debug  
        build/release  
    
and you could configure a debug build in a "debug" subdirectory with

        mkdir build/debug  
        cd build/debug  
        cmake -DCMAKE_BUILD_TYPE=Debug ../../  
        make  
        make install  

To see a detailed compilation log use

        make VERBOSE=1  

To see the available targets use
    
        make help  
    
To build and install the documentation type

        make doc  
        make install  



External libraries
==================

Some functionalities (like solvers ma38 and superlu) require the use of the following libraries:

  - **BLAS**: this library is usually supplied by your computer processor vendor,
            and using a good one is critical to performance.
            If you are unable to locate a vendor BLAS then you should use either
            the GotoBLAS, OpenBLAS or ATLAS BLAS (the latter is often available 
            as part of your linux distribution). If you cannot gain access to any 
            of these, you can obtain the relevant BLAS routines by visiting the 
            following URL: http://www.hsl.rl.ac.uk/blas/
            Such routines obtained from the above url are at least ten times
            slower than the other BLAS libraries mentioned.  
            
  - **LAPACK**: you can obtain the latest LAPACK from http://www.netlib.org/lapack  

  - **SuperLU, SuperLU_MT, SuperLU_DIST**: you can obtain the relevant C functions by visiting the following
            URL: https://portal.nersc.gov/project/sparse/superlu/.
            Note that these three libraries cannot be mixed, only one can be loaded per build target. SuperMU_MT needs either "omp" or "pthread" library. SuperLU_DIST needs "mpi" library.

These librairies are not included in this package. You need to explicitly specify their location in the file "EXTRALIBS". You can also choose to disable the use of these libraries and the associated functionalities by deleting or commenting the locations.


External software
=================

This package contains some fortran routines provided by HSL: "HSL, a collection of Fortran codes for large-scale scientific computation". See http://www.hsl.rl.ac.uk/.


Usage
=====

The bil program prints out the available options when run without any option. To run a specific job, you can enter

        bil [options] myfile  


Bugs/Contact
============

Please mail all bug reports and suggestions to me. I will try to give satisfaction.

email: <patrick.dangla@univ-eiffel.fr>
