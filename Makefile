
all: main.cpp ray.cpp plane.cpp sphere.cpp light.cpp
	g++ -o RayTracer main.cpp ray.cpp plane.cpp sphere.cpp light.cpp -lGL -lGLEW -lGLU -lglut -lX11 -lm -lrt
clean:
	rm RayTracer main.o ray.o sphere.o plane.o
