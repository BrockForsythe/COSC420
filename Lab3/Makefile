main:	main main.o
	mpicc -o main main.o 
main.o: main.c matrix.h
	mpicc -c main.c
clean:
	rm main main.o 
