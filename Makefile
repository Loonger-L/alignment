OMP_PER_READ_THREADS = 1
### use `make OMP_PER_READ_THREADS=N` to override the default 1 OMP thead with N threads

export CPPFLAGS= -g -Wall -O2 -DOMP_PER_READ_THREADS=$(OMP_PER_READ_THREADS) -DHAVE_KALLOC -fopenmp -std=c++11 -Wno-sign-compare -Wno-write-strings -Wno-unused-but-set-variable -fno-tree-vectorize
export LIBS= -lm -lz -lpthread

all:winnowmap

winnowmap: MAKE_DIRS
	+$(MAKE) -e -C src
	$(CXX) $(CPPFLAGS)  src/main.o -o bin/$@ -Lsrc -lwinnowmap $(LIBS)

MAKE_DIRS:
	@if [ ! -e bin ] ; then mkdir -p bin ; fi

clean:
	rm -rf bin
	rm -rf lib
	+$(MAKE) clean -C src

cleanw:
	rm -rf bin/winnowmap
	+$(MAKE) clean -C src
