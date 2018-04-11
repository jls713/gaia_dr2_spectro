## ============================================================================
## General makefile for edf, isochrones and chem
## ============================================================================
## Dependencies:
## 	1. Tact
TACTPATH=/data/jls/tact/tact/

##  2. GSL
GSLPATH =/usr/local/Cellar/gsl/1.16/
ifeq ($(USER),jls)
	GSLPATH =/opt/ioa/software/gsl/2.1/
endif

##  3. GalPot -- implementation in TorusMapper code
TORUSPATH = /home/jls/work/code/Torus_github/Torus/

##  4. nlohmann/json
JSONPATH=/data/jls/json/

##  5. cfitsio, healpix
HEALPIXPATH = /data/jls/Healpix_3.30/

##  6. logging
LOGGPATH = /data/jls/easyloggingpp/

##  7. google test (if want to test)
GTESTPATH = /data/jls/gtest-1.7.0/
####################################################################
## Python module
####################################################################
## These paths are only important if you want to make python modules
## Relies on boost

BOOSTINCPATH =/opt/ioa/software/boost/1.55.0/include
BOOSTLIBPATH =-L/opt/ioa/software/boost/1.55.0/lib
PYTHONINCPATH=/usr/include/python2.7/
ifeq ($(USER),jls)
	PYTHONINCPATH=/opt/ioa/software/python/2.7.8/include/python2.7
endif
PYTHONLIBPATH=-L/opt/ioa/software/python/2.7.8/lib/

## ============================================================================
## ============================================================================
## ============================================================================
CCOMPILER=g++
CCOMPILER_MPI=mpic++
OPTS = --std=c++11 -O3 -fopenmp -fPIC #-ffast-math
# #------------------------------------------------------
# # Set paths to other code
## ============================================================================
## Tact
GENERALPATH = $(TACTPATH)general/
CUBAPATH=/data/jls/tact/tact/general/cuba
CUBAINCPATH=-I$(CUBAPATH)
CUBALIBPATH=-L$(CUBAPATH)
CTPATH = $(GENERALPATH)coordtransforms/
POTPATH = $(TACTPATH)pot/
AAPATH = $(TACTPATH)aa/
GTESTINC=$(GTESTPATH)/include/
GTESTLIB=$(GTESTPATH)/lib/.libs/
## ============================================================================

JSONINC=$(JSONPATH)src/
LOGGINC=$(LOGGPATH)src/

TORUSINCPATH = -I$(TORUSPATH)src/ -I$(TORUSPATH)src/pot/ -I$(TORUSPATH)src/utils/ -I$(TORUSPATH)WDlib/inc/
TORUSLIBPATH = -L$(TORUSPATH)obj/ -L$(TORUSPATH)WDlib/obj/ \
               -L/data/jls/libebf_c_cpp-0.0.3/lib/
HEALPIX_INC = -I$(HEALPIXPATH)include/ \
			  -I$(HEALPIXPATH)src/cxx/generic_gcc/include/ \
			  -I$(HEALPIXPATH)src/cxx/cxxsupport/
HEALPIX_LIB = -L$(HEALPIXPATH)lib/ \
			  -L$(HEALPIXPATH)src/cxx/generic_gcc/lib/ \
			  -Wl,-rpath=$(HEALPIXPATH)src/cxx/generic_gcc/lib/ \
			  -lhealpix_cxx -lchealpix -lcxxsupport -lcfitsio

HDF5PATH=/data/jls/hdf5-1.8.16/
HDF5INCPATH=$(HDF5PATH)include/
HDF5LIBPATH=-L$(HDF5PATH)lib/ -Wl,-rpath=$(HDF5PATH)lib/

INCPY = $(INC) -I$(PYTHONINCPATH) -I$(BOOSTINCPATH)
BOOSTLIB = -lboost_python -lpython2.7
LIBPY = $(PYTHONLIBPATH) $(BOOSTLIBPATH) $(BOOSTLIB)
## ============================================================================
