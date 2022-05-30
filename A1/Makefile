all:
	g++ -fPIC -g -c -Wall psort.cpp -fopenmp
	g++ -shared -Wl,-soname,libsort.so -o libsort.so psort.o -lc
run:
	g++ -fopenmp -Wall psort.cpp driver.cpp -o out
	./out input_1000000.txt 4
clean:
	rm *.so *.o