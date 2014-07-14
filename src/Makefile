ERR = $(shell which clang >/dev/null; echo $$?)
ifneq "$(ERR)" "0"
	CXXFLAGS = -Wall -O3 -fPIC -fforce-addr -funroll-loops -frerun-cse-after-loop -frerun-loop-opt -falign-functions=4
	CPP = g++
else
	CXXFLAGS = -Wall -O3 -dynamiclib -arch x86_64 -funroll-loops -falign-functions=4
	CPP = clang++
endif

OBJECT_FILES = pssm_algorithms.o find_pssm_dna.o

all: find_pssm_dna libpssm.a

find_pssm_dna: libpssm.a find_pssm_dna.o
	$(CPP) find_pssm_dna.o -L. -lpssm -o find_pssm_dna
	
libpssm.a: pssm_algorithms.o
	ar rcs libpssm.a pssm_algorithms.o
	
find_pssm_dna.o: find_pssm_dna.cpp seq_buffer.h
	$(CPP) $(CXXFLAGS) -c find_pssm_dna.cpp -o find_pssm_dna.o

pssm_algorithms.o: pssm_algorithms.cpp pssm_algorithms.hpp seq_buffer.h
	$(CPP) $(CXXFLAGS) -c pssm_algorithms.cpp -o pssm_algorithms.o

clean:
	rm -f  libpssm.a find_pssm_dna $(OBJECT_FILES)
