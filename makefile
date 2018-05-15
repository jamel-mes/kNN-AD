CC=g++
CFLAGS=-Wall -O3
LDLIBS=-lm -lboost_regex

SRC=$(wildcard *.cpp)
OBJ=$(patsubst %.cpp,%.o,$(SRC))
TARGET=kNN

CLEAN=$(OBJ) $(TARGET)


$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) -o  $@ $+ $(LDLIBS)

clean:
	rm $(CLEAN)
