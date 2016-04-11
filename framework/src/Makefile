#Checar en http://www.cs.dartmouth.edu/~campbell/cs50/buildlib.html
#para generacion de una biblioteca estatica a partir de los
#archivos .o generados a partir de la compilacion de un archivos .cpp

#Checar en https://www.gnu.org/software/make/manual/html_node/Wildcard-Function.html
#para analizar el funcionamiento de las funciones "patsubst" y "wildcard"

INCLUDEDIR = .
CC = gcc
CFLAGS=-std=c++0x -I $(INCLUDEDIR)

#modulos a compilar
Containers = $(patsubst containers/%.cpp, ../obj/%.o, $(wildcard containers/*.cpp))
Geometry = $(patsubst geometry/%.cpp, geometry/%.o, $(wildcard geometry/*.cpp))
Image = $(patsubst image/%.cpp, image/%.o, $(wildcard image/*.cpp))
Loaders = $(patsubst loaders/%.cpp, loaders/%.o, $(wildcard loaders/*.cpp))
Matrix = $(patsubst matrix/%.cpp, matrix/%.o, $(wildcard matrix/*.cpp))
Morphing = $(patsubst morphing/%.cpp, morphing/%.o, $(wildcard morphing/*.cpp))
Warping = $(patsubst warping/%.cpp, warping/%.o, $(wildcard warping/*.cpp))

#biblioteca a generar
Lib = morphableKit.a

#objetivo biblioteca
$(Lib) : $(Containers) $(Geometry) $(Image) $(Loaders) $(Matrix) $(Morphing) $(Warping)
	ar  rcs ../obj/lib$(Lib) $(Containers) $(Geometry) $(Image) $(Loaders) $(Matrix) $(Morphing) $(Warping)

#objetivos de los modulos
containers/%.o: containers/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@

geometry/%.o: geometry/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@

image/%.o: image/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@

loaders/%.o: loaders/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@

matrix/%.o: matrix/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@

morphing/%.o: morphing/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@

warping/%.o: warping/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm ../obj/*
