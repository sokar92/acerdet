RELATIVE_PATH= ../..
include $(RELATIVE_PATH)/external.inc

.SUFFIXES:  .o .c .cxx .f .l  .s  .sh  .h  .a

CXX= g++
CLD= $(CXX)
CXXFLAGS = -g -O -ansi -Wall -fPIC -D_REENTRANT $(ROOTINCLUDE)

.cxx.o:
	$(CXX) $(CXXFLAGS) -c $<

all: lib

SLIB1= LibHistoManager
$(SLIB1)_HED= HistoManager.h
$(SLIB1)_OBJ= HistoManager.o	

$(SLIB1)Dict.cxx: $($(SLIB1)_HED) $(SLIB1)LinkDef.h
	@echo "Generating dictionary $@..."
	@rootcint -f $@ -c $^

lib$(SLIB1).so: $(SLIB1)Dict.cxx $($(SLIB1)_OBJ) $(SLIB1)Dict.o
	g++ -shared -o lib$(SLIB1).so -O $($(SLIB1)_OBJ) $(SLIB1)Dict.o

lib$(SLIB1).a: $(SLIB1)Dict.cxx $($(SLIB1)_OBJ) $(SLIB1)Dict.o
	ar cr lib$(SLIB1).a   $($(SLIB1)_OBJ) $(SLIB1)Dict.o

lib: lib$(SLIB1).so lib$(SLIB1).a 

clean:
	rm -rf lib$(SLIB1).so lib$(SLIB1).a $($(SLIB1)_OBJ) $(SLIB1)Dict.*
