#####################################################################
#
#  Name:         Makefile
#  Created by:   Stefan Ritt
#
#  Contents:     Makefile for MIDAS example frontend and analyzer
#
#  $Id$
#
#####################################################################
#
#--------------------------------------------------------------------
# The MIDASSYS should be defined prior the use of this Makefile
ifndef MIDASSYS
missmidas::
	@echo "...";
	@echo "Missing definition of environment variable 'MIDASSYS' !";
	@echo "...";
endif

# get OS type from shell
OSTYPE = $(shell uname)

#--------------------------------------------------------------------
# The following lines contain specific switches for different UNIX
# systems. Find the one which matches your OS and outcomment the 
# lines below.

#-----------------------------------------
# This is for Linux
ifeq ($(OSTYPE),Linux)
OSTYPE = linux
endif

ifeq ($(OSTYPE),linux)
OS_DIR = linux
OSFLAGS = -DOS_LINUX -Dextname
CFLAGS = -g -O2 -Wall
# add to compile in 32-bit mode
# OSFLAGS += -m32
LIBS = -lm -lz -lutil -lnsl -lpthread -lrt
endif

#-----------------------
# MacOSX/Darwin is just a funny Linux
#
ifeq ($(OSTYPE),Darwin)
OSTYPE = darwin
endif

ifeq ($(OSTYPE),darwin)
OS_DIR = darwin
FF = cc
OSFLAGS = -DOS_LINUX -DOS_DARWIN -DHAVE_STRLCPY -DAbsoftUNIXFortran -fPIC -Wno-unused-function
LIBS = -lpthread -lz
LIBS += -lboost_system -lboost_filesystem
SPECIFIC_OS_PRG = $(BIN_DIR)/mlxspeaker
NEED_STRLCPY=
NEED_RANLIB=1
NEED_SHLIB=
NEED_RPATH=

endif

#-----------------------------------------
# ROOT flags and libs
#
ifdef ROOTSYS
ROOTCFLAGS := $(shell  $(ROOTSYS)/bin/root-config --cflags)
ROOTCFLAGS += -DHAVE_ROOT -DUSE_ROOT
ROOTLIBS   := $(shell  $(ROOTSYS)/bin/root-config --libs) -Wl,-rpath,$(ROOTSYS)/lib
ROOTLIBS   += -lThread
else
missroot:
	@echo "...";
	@echo "Missing definition of environment variable 'ROOTSYS' !";
	@echo "...";
endif
#-------------------------------------------------------------------
# The following lines define directories. Adjust if necessary
#                 
DRV_DIR   = $(MIDASSYS)/drivers/camac
INC_DIR   = $(MIDASSYS)/include
LIB_DIR   = $(MIDASSYS)/$(OS_DIR)/lib
SRC_DIR   = $(MIDASSYS)/src

#-------------------------------------------------------------------
# List of analyzer modules
#
MODULES   = movement.o

#-------------------------------------------------------------------
# Hardware driver can be (camacnul, kcs2926, kcs2927, hyt1331)
#
#DRIVER = some driver

#-------------------------------------------------------------------
# Frontend code name defaulted to frontend in this example.
# comment out the line and run your own frontend as follow:
# gmake UFE=my_frontend
#
UFE = ./online/bin/TestMagnetGalilFe
UFE2 = ./online/bin/TestMagnetProbeFe
SRC = ./src/TestMagnetGalilFe
SRC2 = ./src/TestMagnetProbeFe

####################################################################
# Lines below here should not be edited
####################################################################

# MIDAS library
LIB = $(LIB_DIR)/libmidas.a

# compiler
CC = gcc
CXX = g++
CFLAGS += -g -I$(INC_DIR) -I$(DRV_DIR) -I/Users/rhong/Galil/gclib/include
LDFLAGS +=
DYLIBDIR=/Users/rhong/Galil/gclib/dylib

all: $(UFE) $(UFE2)
#all: $(UFE) 

$(UFE): $(LIB) $(LIB_DIR)/mfe.o $(SRC).cpp 
	$(CXX) $(CFLAGS) $(OSFLAGS) -o $(UFE) $(SRC).cpp \
	$(LIB_DIR)/mfe.o $(LIB) \
        $(LDFEFLAGS) $(DYLIBDIR)/gclib.0.dylib $(DYLIBDIR)/gclibo.0.dylib $(LIBS) -std=c++0x
#	$(LDFEFLAGS) $(LIBS) -std=c++0x

$(UFE2): $(LIB) $(LIB_DIR)/mfe.o $(SRC2).cpp 
	$(CXX) $(CFLAGS) $(OSFLAGS) -o $(UFE2) $(SRC2).cpp \
	$(LIB_DIR)/mfe.o $(LIB) \
        $(LDFEFLAGS) $(LIBS) -std=c++0x

clean::
	rm -f *.o *~ \#*
	rm -f $(UFE)
	rm -f $(UFE2)

#end file
