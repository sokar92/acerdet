# -----------------------
# -- compilation flags --
# -----------------------
CC= g++
CFLAGS= -g3 -O2 -fPIC -Wall
HEPMCLIBS= -L$(HEPMCLOCATION)lib -lHepMC

# -------------------------------
# -- Library paths and modules --
# -------------------------------
LIB= AcerDET.so
LIB_MAJOR_VERSION= .2
LIB_MINOR_VERSION= .0

LIB_NAME= lib$(LIB)
LIB_MAJOR_NAME= $(LIB_NAME)$(LIB_MAJOR_VERSION)
LIB_FULL_NAME= $(LIB_MAJOR_NAME)$(LIB_MINOR_VERSION)

INCLUDE_FOLDER= AcerDET
HEADERS_PATH= ./include

PATH_LIB= /usr/lib
PATH_INCLUDE= /usr/include/$(INCLUDE_FOLDER)

SUBTARGETS= module_core \
	module_configuration \
	module_io \
	module_analyse

OBJECTS= AcerDET.o

# ----------------------
# -- Parts of library --
# ----------------------

OBJS= $(shell find . -regex '.*\.o')

# ---------------------------
# -- Instalation procedure --
# ---------------------------

all: prepare compile export_lib export_include
local: prepare compile export_local


# ---------------------------------------------------------------------------------------------------------------------------------------------------
# COMPILATION
# ---------------------------------------------------------------------------------------------------------------------------------------------------

# ---------------------------
# -- Compile whole project --
# ---------------------------

compile: $(SUBTARGETS) $(OBJECTS) $(LIB_FULL_NAME)

# ---------------------------
# -- Create Shared Library --
# -- ------------------------

$(LIB_FULL_NAME): $(OBJS)
	$(CC) -shared $(OBJS) -o $@

# --------------------------
# -- Compile global files --
# --------------------------

%.o : %.cpp
	$(CC) $(CFLAGS) -c $< -o $@ -I$(HEPMCLOCATION)include

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



# ---------------------------------------------------------------------------------------------------------------------------------------------------
# LIBRARY INSTALATION AND FILES MANAGEMENT
# ---------------------------------------------------------------------------------------------------------------------------------------------------

# -----------------
# -- Preparation --
# -----------------

prepare:
	rm -rf $(HEADERS_PATH)
	mkdir $(HEADERS_PATH)
	rm -rf $(INCLUDE_FOLDER)

# --------------------------
# -- Export to global use --
# --------------------------

export_lib:
	cp $(LIB_FULL_NAME) $(PATH_LIB)/$(LIB_FULL_NAME)
	ln -sf $(PATH_LIB)/$(LIB_FULL_NAME) $(PATH_LIB)/$(LIB_MAJOR_NAME)
	ln -sf $(PATH_LIB)/$(LIB_FULL_NAME) $(PATH_LIB)/$(LIB_NAME)
	ldconfig

export_include:
	rm -rf $(PATH_INCLUDE)
	cp -f *.h -t $(HEADERS_PATH)
	cp -r $(HEADERS_PATH) $(INCLUDE_FOLDER)
	mv $(HEADERS_PATH) $(PATH_INCLUDE)

export_local:
	cp -f *.h -t $(HEADERS_PATH)
	mv $(HEADERS_PATH) $(INCLUDE_FOLDER)
	ln -sf $(LIB_FULL_NAME) $(LIB_NAME)

# -------------------
# -- Clean objects --
# -------------------

clean:
	find . -regex '.*\.\(~\|o\|so\)' | xargs rm -rf
	find . -regex '.*\.so\..*' | xargs rm -rf



# ---------------------------------------------------------------------------------------------------------------------------------------------------
# TESTING
# ---------------------------------------------------------------------------------------------------------------------------------------------------

test: prepare compile core_tests


# CORE
core_tests: core_test_histogram

core_test_histogram:
	la
