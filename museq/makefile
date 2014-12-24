###################################################
#
# Makefile for pybam.so
#
###################################################

#
# Macros
#

BOOSTPATH = /share/apps/build/boost_1_51_0
INCLUDE = -I. 
CC = g++
CXX = g++
CXXFLAGS = -O3 -ggdb -DDEBUG_CHECKS $(INCLUDE)
LDFLAGS = -O3 -ggdb 
LDLIBS = 

BINS = \
		pybam.so

SOURCE = \
		src/pybam.cpp

#
# Build
#

all : $(BINS)

pybam.so : $(SOURCE)
		 python setup.py --boost_source=$(BOOSTPATH) install --prefix=./ --install-platlib=./

clean : 
		rm -f *.o *.so

