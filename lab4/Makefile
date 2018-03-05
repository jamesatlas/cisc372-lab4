CC := mpic++
CFLAGS := -O3 -march=native -std=c11

.PHONY: all clean

all: fifteenpuzzle1 fifteenpuzzle2 fifteenpuzzle3

fifteenpuzzle1: fifteenpuzzle1.cpp
	$(CC) -o $@ $^

fifteenpuzzle2: fifteenpuzzle2.cpp
	$(CC) -o $@ $^

fifteenpuzzle3: fifteenpuzzle3.cpp
	$(CC) -o $@ $^

clean:
	rm *.o fifteenpuzzle1 fifteenpuzzle2 fifteenpuzzle3
