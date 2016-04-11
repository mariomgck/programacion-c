LIBDIR = ../../framework/obj
INCLUDEDIR = ../../framework/src/

CC = g++
CFLAGS=-std=c++0x -I $(INCLUDEDIR)
OBJS = $(patsubst %.cpp, %.o, $(wildcard *.cpp))
EXEC = main
LIBS = -L $(LIBDIR) -lmorphableKit -lpthread


$(EXEC) : $(OBJS)
	$(CC) -o $(EXEC) $(OBJS) $(LIBS)

%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm *.o $(EXEC)

run:
	./$(EXEC) ../bootmodel/data/image1.ppm
