

all: RayTracer

RayTracer: main.o ray.o
	g++ -o RayTracer main.o ray.o -lGL -lGLU -lglut -lX11 -lm -lrt

main.o: main.cpp ray.h Angel.h
	g++ -c main.cpp
ray.o: ray.cpp ray.h Angel.h
	g++ -c ray.cpp

clean:
	rm RayTracer main.o ray.o


