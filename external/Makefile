RELATIVE_PATH= ..
include $(RELATIVE_PATH)/common.inc
include $(RELATIVE_PATH)/external.inc

HEPMC_CONVERTER= HepMC_InputConverter
PYTHIA_PROVIDER= Pythia8_ParticleDataProviderFactory
ROOT_HISTOGRAM= Root_HistogramManager
ROOT_NTUPLE= Root_NTupleManager

OBJECTS= $(HEPMC_CONVERTER).o \
	$(PYTHIA_PROVIDER).o \
	$(ROOT_HISTOGRAM).o \
	$(ROOT_NTUPLE).o

LIB= AcerDetExternals
EXT_LIBS= $(PWD)/libHistoManager/libLibHistoManager

# ---------------------
# -- Compile library --
# ---------------------

all: compile_extLibs compile_acerDetLib lib

# --------------------------------
# -- Subtargets ------------------
# --------------------------------
$(HEPMC_CONVERTER).o : $(HEPMC_CONVERTER).cpp
	$(CC) $(CFLAGS) -c $< -o $@ $(HEPMCINCLUDE)

$(PYTHIA_PROVIDER).o : $(PYTHIA_PROVIDER).cpp
	$(CC) $(CFLAGS) -c $< -o $@ $(PYTHIAINCLUDE)

$(ROOT_HISTOGRAM).o : $(ROOT_HISTOGRAM).cpp
	$(CC) $(CFLAGS) -c $< -o $@ $(ROOTINCLUDE)

$(ROOT_NTUPLE).o : $(ROOT_NTUPLE).cpp
	$(CC) $(CFLAGS) -c $< -o $@ $(ROOTINCLUDE)

compile_acerDetLib:
	$(MAKE) -C ../src

compile_extLibs:
	$(MAKE) -C ./libHistoManager

clean_extLibs:
	$(MAKE) -C ./libHistoManager clean

# -------------------------------
# -- Library --------------------
# -------------------------------
lib$(LIB).so: $(OBJECTS)
	$(CC) -shared $(OBJECTS) $(EXT_LIBS).so $(HEPMCLIBS) $(PYTHIALIBS) $(ROOTLIBS) -o lib$(LIB).so

lib: lib$(LIB).so

# ------------------------------
# -- Cleanup -------------------
# ------------------------------
clean: clean_extLibs
	find . -regex '.*\.\(~\|o\|so\)' | xargs rm -rf
