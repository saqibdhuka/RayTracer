#ifndef SPHERE_H
#define SPHERE_H

#include <GL/gl.h>
#include <GL/freeglut.h>
#include <GL/glut.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "Angel.h"
#include "ray.h"

#ifndef MATERIALSTRUCTURE_H
#define MATERIALSTRUCTURE_H

#include "materialStructure.h"

#endif
class Sphere
{
public:

	// int count;
	// struct Material {
	// 	vec3 kdColor;
	// 	vec3 ksColor;
	// 	float specular;
	// 	float reflection;
	// 	float refraction;
	// };

	vec3 center;
	float radius;
	Material mat;


	Sphere();
	Sphere(vec3 point, float rad, Material newMat);
	vec3 getkdColor();
	void draw();
	double intersect(const Ray &ray);
	vec3 get_normal(vec3 rayOrigin);
};


#endif
