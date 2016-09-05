all: makematrix solvematrix calculatebfield
makematrix: makematrix.c
	gcc makematrix.c -o makematrix -lm -fopenmp
solvematrix: solvematrix.o
	gcc -fopenmp -o solvematrix solvematrix.o -L~/lis-1.5.65/lib -llis -lm
solvematrix.o: solvematrix.c
	gcc -lm -c -fopenmp -I~/lis-1.5.65/include solvematrix.c
calculatebfield: calculatebfield.cpp
	g++ -o calculatebfield calculatebfield.cpp -pthread -lm `root-config --cflags --libs`
