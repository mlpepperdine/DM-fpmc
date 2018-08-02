####################################################
# M. Saimpert-interface for dm_dm coupling 
####################################################

.PHONY: clear test
.SUFFIXES: .cpp .o

CC=g++
CFLAGS=-g -Wall -DNUM_DOUBLE 
SOFLAGS=-fPIC
GSL=`gsl-config --cflags --libs`

all: dm_sqme_dm.a libdm_sqme_dm.so

objects=functions.o helicity_amplitudes.o sqme.o
f_0_objects=$(objects:%=src/%)

dm_sqme_dm.a:$(f_0_objects)
	ar -r $@ $(f_0_objects)

libdm_sqme_dm.so:$(f_0_objects)
	$(CC) $(SOFLAGS) -shared $(f_0_objects) -o $@ $(GSL)

.cpp.o:
	$(CC) $(CFLAGS) $(SOFLAGS) $(INS) -c $< -o $@ $(GSL)

clean:
	rm -rf src/*.o test/*.o *.a *.so testcall testcall2

