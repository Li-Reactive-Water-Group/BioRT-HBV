#-----------------------------------------------------------------
# BioRT Makefile
# -----------------------------------------------------------------

CC = gcc
CFLAGS = -g -O0 -Wall -Wextra

CMAKE_VER_NUM := $(shell cmake --version 2> /dev/null |awk '{print $$3}')
ifeq ($(CMAKE_VER_NUM),)
	CMAKE_VER_NUM := 0.0.0
endif
CMAKE_REQ_VER = 3.1.3
CMAKETEST := $(shell printf '%s\n' $(CMAKE_VER_NUM) $(CMAKE_REQ_VER) | sort -V | head -n 1)

ifeq ($(CMAKETEST),$(CMAKE_REQ_VER))
	CMAKE_EXIST = 1
	CMAKE=cmake
else
	CMAKE_EXIST = 0
	OS := $(shell uname)
	ifeq ($(OS),Darwin)
		CMAKE_VERS = cmake-3.7.2-Darwin-x86_64
		CMAKE = $(PWD)/$(CMAKE_VERS)/CMake.app/Contents/bin/cmake
	else
		CMAKE_VERS = cmake-3.7.2-Linux-x86_64
		CMAKE = $(PWD)/$(CMAKE_VERS)/bin/cmake
	endif
endif

CVODE_PATH = ./cvode/instdir
CVODE_LIB = $(CVODE_PATH)/lib

SRCDIR = ./src
LIBS = -lm
INCLUDES = \
	-I$(SRCDIR)/include\
	-I$(CVODE_PATH)/include

LFLAGS = -lsundials_cvode -L$(CVODE_LIB) -lsundials_nvecserial

SRCS_ = main.c\
	custom_io.c\
	init.c\
	lookup.c\
	optparse.c\
	print.c\
	react.c\
	read_chem.c\
	read_cini.c\
  read_precipchem.c\
	read_hbv.c\
	read_param.c\
	read_soil.c\
	set_numexp.c\
	speciation.c\
	time_func.c\
	transpt.c\
	util_func.c

HEADERS_ = include/biort.h\
	include/custom_io.h

EXECUTABLE = biort
MSG = "...  Compiling BioRT  ..."

SRCS = $(patsubst %,$(SRCDIR)/%,$(SRCS_))
HEADERS = $(patsubst %,$(SRCDIR)/%,$(HEADERS_))
OBJS = $(SRCS:.c=.o)

.PHONY: all clean cvode cmake

biort: $(OBJS)
	@echo
	@echo $(MSG)
	@echo
	@$(CC) $(CFLAGS) $(SFLAGS) $(INCLUDES) -o $(EXECUTABLE) $(OBJS) $(LFLAGS) $(LIBS)

all:    cvode biort

cmake:
ifneq ($(CMAKE_EXIST),1)
	@echo "CVODE installation requires CMake v$(CMAKE_REQ_VER) or above."
	@echo "Download CMake $(CMAKE_VERS) from cmake.org"
	@curl https://cmake.org/files/v3.7/$(CMAKE_VERS).tar.gz -o $(CMAKE_VERS).tar.gz &> /dev/null
	@echo
	@echo "Extract $(CMAKE_VERS).tar.gz"
	@tar xzf $(CMAKE_VERS).tar.gz
endif

cvode:	cmake
	@echo "Install CVODE library"
	@cd cvode && mkdir -p instdir && mkdir -p builddir
	@cd $(CVODE_PATH) && $(CMAKE) -DCMAKE_INSTALL_PREFIX=../instdir -DCMAKE_INSTALL_LIBDIR=lib -DBUILD_SHARED_LIBS=OFF -DEXAMPLES_ENABLE_C=OFF -DEXAMPLES_INSTALL=OFF ../
	@cd $(CVODE_PATH) && make && make install
	@echo "CVODE library installed."
ifneq ($(CMAKE_EXIST),1)
	@echo "Remove CMake files"
	@$(RM) -r $(CMAKE_VERS).tar.gz $(CMAKE_VERS)
endif

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) $(SFLAGS) $(INCLUDES) -c $<  -o $@


clean:			## Clean executables and objects
	@echo
	@echo "... Cleaning ..."
	@echo
	@$(RM) $(SRCDIR)/*.o *~ $(EXICUTABLE)
