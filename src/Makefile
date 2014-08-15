LIBDIR= $(shell root-config --libdir)
INCDIR= $(shell root-config --incdir)
#  Root  libraries for linker
RLIB=   -L$(LIBDIR) -lCore -lCint -lHist -lGraf -lGraf3d -lGpad -lTree -lRint \
	-lPostscript -lMatrix -lPhysics -lm -ldl -lpthread -rdynamic
# -----------------------
# -- compilation flags --
# -----------------------
CC= g++
CFLAGS= -g3 -O2 -fPIC -Wall
HEPMCLIBS= -L$(HEPMCLOCATION)lib -lHepMC
ROOTLIBS := $(shell root-config --libs)
PWD := $(shell pwd)
HISTOLOCATION= $(PWD)/../libHistoManager

# -------------------------------
# -- Library paths and modules --
# -------------------------------
LIB= AcerDET
LIB_MAJOR_VERSION= .2
LIB_MINOR_VERSION= .0

LIB_NAME= lib$(LIB).so
LIB_MAJOR_NAME= $(LIB_NAME)$(LIB_MAJOR_VERSION)
LIB_FULL_NAME= $(LIB_MAJOR_NAME)$(LIB_MINOR_VERSION)

HEADERS_PATH= ./include

SUBTARGETS= module_core \
	module_configuration \
	module_io \
	module_analyse \
	module_external

OBJECTS= AcerDET.o

# ----------------------
# -- Parts of library --
# ----------------------

OBJS= $(shell find . -regex '.*\.o')

# ---------------------------
# -- Instalation procedure --
# ---------------------------

all: prepare compile copy_headers link_lib

# ---------------------------------------------------------------------------------------------------------------------------------------------------
# COMPILATION
# ---------------------------------------------------------------------------------------------------------------------------------------------------

# ---------------------------
# -- Compile whole project --
# ---------------------------

compile: $(SUBTARGETS) $(OBJECTS) $(LIB_FULL_NAME)

# ---------------------------
# -- Create Shared Library --
# ---------------------------

$(LIB_FULL_NAME): $(OBJS)
	$(CC) -shared $(HISTOLOCATION)/libLibHistoManager.so $(OBJS) -o $@ $(ROOTLIBS)

# --------------------------
# -- Compile global files --
# --------------------------

%.o : %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

# -------------------------
# -- Library sub modules --
# -------------------------

module_core:
	$(MAKE) -C ./core

module_configuration:
	$(MAKE) -C ./conf

module_io: 
	$(MAKE) -C ./io
	
module_analyse:
	$(MAKE) -C ./analyse

module_external:
	$(MAKE) -C ./external

# ---------------------------------------------------------------------------------------------------------------------------------------------------
# LIBRARY INSTALATION AND FILES MANAGEMENT
# ---------------------------------------------------------------------------------------------------------------------------------------------------

prepare:
	rm -rf $(HEADERS_PATH)
	mkdir $(HEADERS_PATH)

link_lib:
	ln -sf $(LIB_FULL_NAME) $(LIB_MAJOR_NAME)
	ln -sf $(LIB_FULL_NAME) $(LIB_NAME)

copy_headers:
	cp -f *.h -t $(HEADERS_PATH)

# -------------------
# -- Clean objects --
# -------------------

clean:
	find . -regex '.*\.\(~\a|\|o\|so\)' | xargs rm -rf
	find . -regex '.*\.so\..*' | xargs rm -rf
	rm -rf $(HEADERS_PATH)
