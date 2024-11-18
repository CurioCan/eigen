#!/bin/bash

# Specify the compiler (adjust as needed)
CC=gcc

# Compile the C code
$CC -Wall -Werror -lm -o code code.c -lm

# Check if compilation was successful
if [ $? -ne 0 ]; then
  echo "Compilation failed. Please check for errors."
  exit 1
fi

# Run the executable
./code

# Optionally, remove the executable
rm code
