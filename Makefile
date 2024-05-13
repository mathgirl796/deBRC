
SRC := $(wildcard src/*.cpp) $(wildcard src/kmc_api/*.cpp) $(wildcard src/klib/*.cpp) $(wildcard src/FastaReader/*.cpp)
LIB := -Wl,-Bstatic -Llib -Wl,-Bdynamic -lkmc_core -lbz2 -lz -lpthread
OBJ := $(patsubst src/%.cpp, obj/%.o, $(SRC))
OUT := bin

GCC := g++
DEBUG := -DDEBUG 
CFLAGS := -Wall -g -Wl,--allow-multiple-definition 

TARGET := $(OUT)/main

$(TARGET): $(OBJ)
	@mkdir -p $(shell dirname $@)
	$(GCC) -o $@ $(OBJ) $(LIB) $(DEBUG) $(CFLAGS)

obj/%.o: src/%.cpp
	@mkdir -p $(shell dirname $@)
	$(GCC) -o $@ -c $< $(DEBUG) $(CFLAGS)

clean:
	rm -rf bin/*
	rm -rf obj/*

clean_main:
	rm $(TARGET)

rebuild: clean $(TARGET)
main: clean_main $(TARGET)