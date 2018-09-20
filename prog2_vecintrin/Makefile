all: myexp

logger.o: logger.cpp logger.h CS348Kintrin.h CS348Kintrin.cpp
	g++ -c logger.cpp

CS348Vintrin.o: CS348Kintrin.cpp CS348Kintrin.h logger.cpp logger.h
	g++ -c CS348Kintrin.cpp

myexp: CS348Kintrin.o logger.o main.cpp
	g++ -I../common logger.o CS348Kintrin.o main.cpp -o myexp

clean:
	rm -f *.o myexp
