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

PATH_LIB= /usr/lib
PATH_INCLUDE= /usr/include/AcerDET

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
	rm -rf ./include
	mkdir include

# --------------------------
# -- Export to global use --
# --------------------------

export_lib:
	cp $(LIB_FULL_NAME) $(PATH_LIB)/$(LIB_FULL_NAME)
	ln -sf $(PATH_LIB)/$(LIB_FULL_NAME) $(PATH_LIB)/$(LIB_MAJOR_NAME)
	ln -sf $(PATH_LIB)/$(LIB_FULL_NAME) $(PATH_LIB)/$(LIB_NAME)
	ldconfig

HEADERS_PATH= ./include
export_include:
	rm -rf $(PATH_INCLUDE)
	cp -f *.h -t $(HEADERS_PATH)
	mv ./include  $(PATH_INCLUDE)

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
