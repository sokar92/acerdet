RELATIVE_PATH= .
include common.inc
include external.inc

DEMO= demo.exe demoCreateConfig.exe

EXT_LIB= $(PWD)/$(RELATIVE_PATH)/external/libAcerDetExternals.so
ACERDET_LIB= $(PWD)/$(RELATIVE_PATH)/src/libAcerDET.so

all: acerDetDeps demo

acerDetDeps:
	$(MAKE) -C ./external

demo:
	$(CC) $(CFLAGS) -c demo.cpp -o demo.o $(PYTHIAINCLUDE) $(HEPMCINCLUDE) -D_REENTRANT $(ROOTINCLUDE) -I ./src
	$(CC) $(CFLAGS) -c demoCreateConfig.cpp -o demoCreateConfig.o $(PYTHIAINCLUDE) $(HEPMCINCLUDE) -I ./src
	$(CC) demo.o -o demo.exe $(ACERDET_LIB) $(EXT_LIB) $(PYTHIALIBS) $(HEPMCLIBS) $(ROOTLIBS) -Wl,-rpath,$(PYTHIALOCATION)lib -Wl,-rpath,$(HEPMCLOCATION)lib
	$(CC) demoCreateConfig.o -o demoCreateConfig.exe $(ACERDET_LIB) $(EXT_LIB) $(PYTHIALIBS) $(HEPMCLIBS) -Wl,-rpath,$(PYTHIALOCATION)lib -Wl,-rpath,$(HEPMCLOCATION)lib

clean:
	rm -f *~ *.o *.exe *.a
