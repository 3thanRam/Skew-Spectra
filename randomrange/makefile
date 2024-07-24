CC = gcc
CFLAGS  = -Wall 
LIBFLAGS= -lm -ldl

OBJ_FILES=output/Main.o output/randomrange.o output/ini.o

all : Main
Main : $(OBJ_FILES)
	$(CC) $(CFLAGS) -o $@ $^  $(LIBFLAGS)
output/%.o: src/%.c
	$(CC) $(CFLAGS) -c $< -o $@
clean:
	rm -f output/*.o Main


