CC = gcc
CFLAGS = -Iinclude -Wall -Wextra -Werror -fPIC -O2

SRC_DIR = src
SRC = $(wildcard $(SRC_DIR)/*.c)
# SRC = src/mat_utils.c src/matrix.c
OBJ = $(SRC:.c=.o)

TEST_DIR = test
TEST_SRC = $(wildcard $(TEST_DIR)/*.c)
TEST_OBJ = $(TEST_SRC:.c=.o)
TEST_TARGETS = $(TEST_OBJ:.o=)

MAIN = main.c
MAIN_OBJ = $(MAIN:.c=.o)
TARGET = main

# LIB = libQPIPM.so
LIB_DIR = .
# LIBS = -lmatrix
LIBS = -lm

all: $(TARGET) $(TEST_TARGETS)

# Compile source files to object files
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# # Build the shared library
# $(LIB): $(OBJ)
# 	$(CC) -shared -o $(LIB) $(OBJ)

# Link the main program
$(TARGET): $(OBJ) $(MAIN_OBJ)
	$(CC) $(OBJ) $(MAIN_OBJ) -L$(LIB_DIR) $(LIBS) -o $@

# Compile and link test programs
$(TEST_DIR)/%: $(TEST_DIR)/%.o $(OBJ)
	$(CC) $(OBJ) $< -L$(LIB_DIR) $(LIBS) -o $@

clean:
	rm -f $(OBJ) $(MAIN_OBJ) $(TEST_OBJ) $(TARGET) $(TEST_TARGETS)

.PHONY: all clean