# --------------------------------------------
# -- User changeable region ------------------
# --------------------------------------------

HEPMCLOCATION= ../pythia8-hepmc-example/HepMC-2.06.09/
PYTHIALOCATION= ../pythia8-hepmc-example/pythia8185/
ROOTLOCATION= ../root/

# --------------------------------------------
# -- DO NOT CHANGER ANYTHING BELOW -----------
# --------------------------------------------
PWD= $(shell pwd)

HEPMCINCLUDE= -I$(PWD)/$(RELATIVE_PATH)/$(HEPMCLOCATION)include
PYTHIAINCLUDE= -I$(PWD)/$(RELATIVE_PATH)/$(PYTHIALOCATION)include
ROOTINCLUDE= -I$(PWD)/$(RELATIVE_PATH)/$(ROOTLOCATION)include

HEPMCLIBS= -L$(PWD)/$(RELATIVE_PATH)/$(HEPMCLOCATION)lib -lHepMC
PYTHIALIBS= -L$(PWD)/$(RELATIVE_PATH)/$(PYTHIALOCATION)lib -lpythia8 -llhapdfdummy -lpythia8tohepmc

ROOTLIBS := $(shell root-config --libs) -lCintex
