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
		if [ -z "$(BOOSTPATH)" ] || [ -z "$(PYTHON)" ] ; then echo "please set the boost and python paths i.e. run make PYTHON=/path/to/python BOOSTPATH=/path/to/boost"; exit -1; fi
		$(PYTHON) setup.py --boost_source=$(BOOSTPATH) install --prefix=./ --install-platlib=./

clean :
		rm -rf *.o *.so build
