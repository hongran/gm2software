all: preparemakematrix makematrix solvematrix bfield bfieldHong
preparemakematrix: preparemakematrix.c
	gcc preparemakematrix.c -o preparemakematrix -lm -fopenmp
makematrix: makematrix.c
	gcc makematrix.c -o makematrix -lm -fopenmp
solvematrix: solvematrix.o
	gcc -fopenmp -o solvematrix solvematrix.o -L~/lis-1.5.65/lib -llis -lm
solvematrix.o: solvematrix.c
	gcc -lm -c -fopenmp -I~/lis-1.5.65/include solvematrix.c
bfield: bfield.c
	gcc -o bfield bfield.c -lm
bfieldHong: bfieldHong.c
	gcc -o bfieldHong bfieldHong.c -lm
