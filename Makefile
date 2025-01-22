# Variables
CC = mpicc
CFLAGS = -Wall -O2
LDFLAGS = -lm
TARGET = bin/mandel
SRC_DIR = src
SRC = $(wildcard $(SRC_DIR)/*.c)
OBJ = $(SRC:.c=.o)

# Default target
all: $(TARGET)

# Compile the program
$(TARGET): $(OBJ)
	mkdir -p bin
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# Clean up
clean:
	rm -f $(SRC_DIR)/*.o $(TARGET)

# Phony targets
.PHONY: all clean

