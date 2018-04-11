## ============================================================================
## Makefile for edf code
## ============================================================================
include ../Makefile.inc
# #------------------------------------------------------
# # Set include and library paths
ISOPATH=/home/jls/work/code/gaia_sprint/isochrone/
EDFPATH=/home/jls/work/code/gaia_sprint/edf/
INC = -I../ -Iinc/ -I$(GENERALPATH) -I$(CTPATH)inc/ \
	  -I$(POTPATH)inc/ -I$(AAPATH)inc/ $(TORUSINCPATH) -I$(JSONINC) \
	  -I$(GSLPATH)include/ -I$(CTPATH)inc/ -I$(HEALPIX_INC) \
	   -I$(ISOPATH)inc/  -I$(EDFPATH)inc/ -I$(HDF5INCPATH) -I$(LOGGINC) \
	   -I$(CUBAPATH) -I../chemevo/inc/ -I$(GTESTINC)

LIBPATH = -L$(EDFPATH)lib/ -L$(ISOPATH)lib/ $(TORUSLIBPATH) -L$(AAPATH)lib/ \
		  -L$(POTPATH) -L$(CUBAPATH) -L$(GSLPATH)lib/ -L$(CTPATH) \
		  $(HEALPIX_LIB) -Wl,-rpath=$(EDFPATH)lib/ -Wl,-rpath=$(ISOPATH)lib/ \
		  $(HDF5LIBPATH)
LIBS = -lisodist_js -ledf -lcuba -lcoords -laa -lpot_js -lTorus -lPot -lOther \
	   -lWD -lebf_cpp -lhdf5_cpp -lhdf5 -lgsl -lgslcblas

LIB = $(LIBPATH) $(LIBS)
#------------------------------------------------------
# Set files
SRC = src/
MAINSRC = src/mains/
MAINEXE = mains/
OBJ = obj/
LIBDIR = lib/
PY = py/

OBJECTS = $(OBJ)sf.o $(OBJ)sf_impl.o $(OBJ)mass_function.o \
		  $(OBJ)ages.o $(OBJ)imf.o $(OBJ)params.o

#------------------------------------------------------

default: all

all: $(PY)edf_sampling.so

forceall:	clean all

$(OBJ)%.o : ../chemevo/src/%.cpp
	$(CCOMPILER) $(OPTS) -c $< -o $@ $(INC) -DTORUS -DGSL_2

$(OBJ)%.o : $(SRC)%.cpp
	$(CCOMPILER) $(OPTS) -c $< -o $@ $(INC) -DTORUS -DGSL_2

$(OBJ)%.o : $(MAINSRC)%.cpp
	$(CCOMPILER) $(OPTS) -c $< -o $@ $(INC) -DTORUS -DGSL_2

$(MAINEXE)%.exe : $(OBJ)%.o $(OBJECTS)
	$(CCOMPILER) $(OPTS) $(OBJECTS) $< -o $@ $(LIB) -DTORUS -DGSL_2

test.exe: $(OBJECTS) obj/test_edf.o
	$(CCOMPILER) $(OPTS) $(OBJECTS) obj/test_edf.o $(LIB) -L$(GTESTLIB) -Wl,-R$(GTESTLIB) -lgtest -o test.exe -DGSL_2 -DTORUS
clean:
	rm -rf obj/*.o; rm -rf *.so; rm -rf *.exe $(ALL)

####################################################################
$(OBJ)/edf_sampling.o : $(SRC)/edf_sampling.cpp
	$(CCOMPILER) --std=c++11 $(OPTS) -c $< -o $@ $(INCPY) -DTORUS -DGSL_2

$(PY)edf_sampling.so : $(OBJ)/edf_sampling.o $(OBJECTS)
	$(CCOMPILER) -shared -Wl,--no-undefined -Wl,-z,defs $(OPTS) $(OBJECTS) $(OBJ)/edf_sampling.o $(LIBPY) $(LIB) -o $(PY)edf_sampling.so -DTORUS
