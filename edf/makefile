## ============================================================================
## Makefile for edf code
## ============================================================================
include ../Makefile.inc
## ============================================================================
# # Set include and library paths
INC = -Iinc/ -I$(GENERALPATH) -I$(CTPATH)inc/ \
	  -I$(POTPATH)inc/ -I$(AAPATH)inc/ $(TORUSINCPATH) -I$(JSONINC) \
	  -I$(GSLPATH)include/ -I$(CTPATH)inc/
LIBPATH = $(TORUSLIBPATH) -L$(AAPATH)lib/ -L$(POTPATH) -L$(CUBAPATH) \
		  -L$(GSLPATH)lib/ -L$(CTPATH)
LIBS = -lcuba -lcoords -laa -lpot_js -lTorus -lPot -lOther -lWD -lebf_cpp \
       -lgsl -lgslcblas
LIB = $(LIBS) $(LIBPATH)
#------------------------------------------------------
# Set files
LIBDIR = lib/
SRC = src/
MAINSRC = src/mains/
MAINEXE = mains/
OBJ = obj/
PY = py/

OBJECTS = $(OBJ)sfr.o $(OBJ)metal_relation.o $(OBJ)thin_disc.o \
		  $(OBJ)thick_disc.o $(OBJ)halo.o
MAINOBJECTS = $(OBJECTS) $(OBJ)edf.o

#------------------------------------------------------

default: all

all:  mains lib

forceall:	clean all

mains: $(MAINEXE)plotFeH.exe $(MAINEXE)plotGR.exe $(MAINEXE)plotRotCurve.exe \
	   $(MAINEXE)plotUVW.exe $(MAINEXE)global_metalsplit.exe \
	   $(MAINEXE)global.exe $(MAINEXE)global_1D.exe $(LIBDIR)libedf.so

lib: $(LIBDIR)libedf.so

$(OBJ)%.o : $(SRC)%.cpp
	$(CCOMPILER) $(OPTS) -c $< -o $@ $(INC) -DTORUS

$(OBJ)%.o : $(MAINSRC)%.cpp
	$(CCOMPILER) $(OPTS) -c $< -o $@ $(INC) -DTORUS

$(MAINEXE)%.exe : $(OBJ)%.o $(MAINOBJECTS)
	$(CCOMPILER) $(OPTS) $(MAINOBJECTS) $< -o $@ $(LIB) -DTORUS

$(LIBDIR)libedf.so : $(OBJ)edf.o $(OBJECTS)
	$(CCOMPILER) -shared -Wl,--no-undefined -Wl,-z,defs $(OPTS) $(OBJECTS) $< -o $@ $(LIB) -DTORUS

clean:
	rm -rf obj/*.o; rm -rf *.so; rm -rf *.exe $(ALL)

####################################################################

$(OBJ)/edf_py.o : $(PY)/edf_py.cpp
	$(CCOMPILER) $(OPTS) -c $< -o $@ $(INCPY) -DTORUS

$(PY)edf_py.so : $(OBJ)/edf_py.o $(MAINOBJECTS)
	$(CCOMPILER) -shared -Wl,--no-undefined -Wl,-z,defs $(OPTS) $(MAINOBJECTS) $(OBJ)/edf_py.o -o $(PY)edf_py.so $(LIBPY) $(LIB) -DTORUS
