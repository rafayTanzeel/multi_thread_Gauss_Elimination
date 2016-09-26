all: gauss

CC = gcc
CCFLAGS = -g -Wall -O3 -std=c11
LDFLAGS = -g -Wall -O3
LFLAGS = -lm -lpthread

gauss: gauss.o
	$(CC) $(LDFLAGS) -o gauss gauss.o $(LFLAGS)

gauss.o: gauss.c
	$(CC) $(CCFLAGS) -c gauss.c -o gauss.o

clean:
	-rm *.o *~ gauss
