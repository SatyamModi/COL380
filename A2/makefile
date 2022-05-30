compile:
	mpic++ -std=c++17 -o wtf *.cpp librandomizer.a
run:
	mpirun -np 1 ./wtf data/81867_545671/edges.dat 81867 545671 0.1 250 30 20 369
clean:
	rm wtf