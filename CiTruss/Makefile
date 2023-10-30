UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	# Mac OSX
	METISFILE=metis
	CXX=clang++
	CXXFLAGS=-fopenmp -fPIC -O3 -DNDEBUG
	INCLUDEFLAGS=-I/usr/local/opt/llvm/include -I../metis-5.1.0/include -std=c++11
	LIBFLAGS=-L/usr/local/opt/llvm/lib -L../metis-5.1.0/lib -lm
else
	# GNU/Linux platform
	METISFILE=altmetis
	CXX=g++
	CXXFLAGS=-fopenmp -fPIC -O3 -DNDEBUG
	INCLUDEFLAGS=-I/usr/include -I/usr/local/include -I../metis-5.1.0/include -std=c++11
	LIBFLAGS=-L/usr/local/lib -L../metis-5.1.0/lib -lm
endif



all: citruss_sum citruss_diff mega_scggm.o sparse.o smat.o util.o

citruss_sum: citruss_sum_run.cpp mega_scggm.o smat.o sparse.o util.o
	$(CXX) ${CXXFLAGS} ${INCLUDEFLAGS} ${LIBFLAGS} $^ -l${METISFILE} -o citruss_sum

citruss_diff: citruss_diff_run.cpp mega_scggm.o smat.o sparse.o util.o
	$(CXX) ${CXXFLAGS} ${INCLUDEFLAGS} ${LIBFLAGS} $^ -l${METISFILE} -o citruss_diff

mega_scggm.o: mega_scggm.h mega_scggm.cpp sparse.h smat.h util.h
	${CXX} -c ${CXXFLAGS} ${INCLUDEFLAGS} ${LIBFLAGS} mega_scggm.cpp -o mega_scggm.o

sparse.o: sparse.h sparse.cpp
	${CXX} -c ${CXXFLAGS} ${INCLUDEFLAGS} ${LIBFLAGS} sparse.cpp -o sparse.o

smat.o: smat.h smat.cpp
	${CXX} -c ${CXXFLAGS} ${INCLUDEFLAGS} ${LIBFLAGS} smat.cpp -o smat.o

util.o: util.h util.cpp
	${CXX} -c ${CXXFLAGS} ${INCLUDEFLAGS} ${LIBFLAGS} util.cpp -o util.o

clean:
	rm -rf *.pyc *.o citruss_sum citruss_diff *.mexa64 *.mexmaci64
