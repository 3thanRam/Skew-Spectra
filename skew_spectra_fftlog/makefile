CC = gcc
CFLAGS  = -Wall 
LIBFLAGS= -lgsl -lfftw3 -lm -ldl -lblas  -fPIC

OBJ_FILES=output/Main.o output/Util.o output/Plotting.o output/fftw_funcs.o output/class_funcs.o output/ini.o

all : Main
Main : $(OBJ_FILES)
	$(CC) $(CFLAGS) -o $@ $^  $(LIBFLAGS)
output/%.o: src/%.c
	$(CC) $(CFLAGS) -c $< -o $@
clean:
	rm -f output/*.o Main