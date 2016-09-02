RELATIVE_PATH= .
include common.inc
include external.inc

DEMO= demo.exe demoCreateConfig.exe

EXT_LIB= $(PWD)/$(RELATIVE_PATH)/external/libAcerDetExternals.so
ACERDET_LIB= $(PWD)/$(RELATIVE_PATH)/src/libAcerDET.so

EVENTS= 100000

T1= hadronisation
T2= hardprocess
T3= IsrFsr

DIR= res
DIR1= $(DIR)/$(T1)
DIR2= $(DIR)/$(T2)
DIR3= $(DIR)/$(T3)

all: acerDetDeps demo

acerDetDeps:
	$(MAKE) -C ./external

demo:
	$(CC) $(CFLAGS) -c demo.cpp -o demo.o $(PYTHIAINCLUDE) $(HEPMCINCLUDE) -D_REENTRANT $(ROOTINCLUDE) -I ./src
	$(CC) $(CFLAGS) -c demoCreateConfig.cpp -o demoCreateConfig.o $(PYTHIAINCLUDE) $(HEPMCINCLUDE) -I ./src
	$(CC) demo.o -o demo.exe $(ACERDET_LIB) $(EXT_LIB) $(PYTHIALIBS) $(HEPMCLIBS) $(ROOTLIBS) -Wl,-rpath,$(PYTHIALOCATION)lib -Wl,-rpath,$(HEPMCLOCATION)lib $(FASTJETLIBS)
	$(CC) demoCreateConfig.o -o demoCreateConfig.exe $(ACERDET_LIB) $(EXT_LIB) $(PYTHIALIBS) $(HEPMCLIBS) -Wl,-rpath,$(PYTHIALOCATION)lib -Wl,-rpath,$(HEPMCLOCATION)lib $(FASTJETLIBS)

clean:
	rm -f *~ *.o *.exe *.a

run_cleanUp:
	rm -rf $(DIR)

run_prepareDirectory: run_cleanUp
	mkdir $(DIR)
	mkdir $(DIR1)
	mkdir $(DIR2)
	mkdir $(DIR3)

run: demo.exe run_prepareDirectory
	cp acerdet_old.dat acerdet.dat
	./demo.exe conf/$(T1).conf $(EVENTS)
	cp conf/$(T1).conf.root $(DIR1)/old.root
	cp acerdet_fast.dat acerdet.dat
	./demo.exe conf/$(T1).conf $(EVENTS)
	cp conf/$(T1).conf.root $(DIR1)/fast.root
	cp acerdet_old.dat acerdet.dat
	./demo.exe conf/$(T2).conf $(EVENTS)
	cp conf/$(T2).conf.root $(DIR2)/old.root
	cp acerdet_fast.dat acerdet.dat
	./demo.exe conf/$(T2).conf $(EVENTS)
	cp conf/$(T2).conf.root $(DIR2)/fast.root
	cp acerdet_old.dat acerdet.dat
	./demo.exe conf/$(T3).conf $(EVENTS)
	cp conf/$(T3).conf.root $(DIR3)/old.root
	cp acerdet_fast.dat acerdet.dat
	./demo.exe conf/$(T3).conf $(EVENTS)
	cp conf/$(T3).conf.root $(DIR3)/fast.root
