IFLAGS=-I/opt/homebrew/include
LFLAGS=-L/opt/homebrew/lib -L/opt/homebrew/opt/openblas/lib -lm -lgsl -lgslcblas
FLAGS=-g -pg -Wall -Wextra -Wpedantic --std=c99
CC=gcc $(FLAGS) $(IFLAGS)

main: main.c bon16.o cap09.o ued14.o laf05.o
	$(CC) $^ -o $@ $(LFLAGS)

dz.h: create_dz.py
	python3 create_dz.py

%.o: %.c %.h
	$(CC) -c $<


test_cap09: test_cap09.o cap09.o
	$(CC) -o $@ $(LFLAGS)

test_laf05: test_laf05.c laf05.o
	$(CC) $^ -o $@ $(LFLAGS)
