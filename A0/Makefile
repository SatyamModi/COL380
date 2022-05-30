CC=g++
#g++ -std=c++11 -O2 -mavx2 -c classify.cpp 
CFLAGS=-std=c++11 -O2 -g 

sources=main.cpp classify.cpp
objects=$(sources:.cpp=.o)

classify:$(objects)
	$(CC) $(CFLAGS) -fopenmp $^ -o $@

run: classify
	./classify rfile dfile 1009072 4 3

#classify.o: classify.h
#	touch classify.o
#	$(CC) $(CFLAGS) -c $<
#main.o: classify.h
#	$(CC) $(CFLAGS) -c $<

%.o: %.cpp classify.h
	$(CC) $(CFLAGS) -fopenmp -c $<
