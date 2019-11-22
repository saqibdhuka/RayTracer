#include "sphere.h"
#include <GL/glew.h>
#include <GL/gl.h>
#include <GL/freeglut.h>
#include <GL/glut.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "Angel.h"




Sphere :: Sphere(){
	center.x = 0.0;
	center.y = 0.0;
	center.z = 0.0;

	radius = 0.0;
	
	mat.kdColor.x = 0.0;
	mat.kdColor.y = 0.0;
	mat.kdColor.z = 0.0;

	mat.ksColor.x = 0.0;
	mat.ksColor.y = 0.0;
	mat.ksColor.z = 0.0;

	mat.specular = 0.0;
	mat.reflection = 0.0;
	mat.refraction = 0.0;

	
}

Sphere:: Sphere(vec3 point, float rad, Material newMat){
	center.x = point.x;
	center.y = point.y;
	center.z = point.z;

	radius = rad;
	
	newMat.kdColor.x = mat.kdColor.x;
	newMat.kdColor.y = mat.kdColor.y;
	newMat.kdColor.z = mat.kdColor.z;

	newMat.ksColor.x = mat.ksColor.x;
	newMat.ksColor.y = mat.ksColor.y;
	newMat.ksColor.z = mat.ksColor.z;

	newMat.specular = mat.specular;
	newMat.reflection = mat.reflection;
	newMat.refraction = mat.refraction;
}

void Sphere :: draw(){
	glutSolidSphere(radius, 20.0, 50.0);
}


vec3 Sphere::getkdColor(){
	return mat.kdColor;
}

bool Sphere :: insersect(const Ray &ray, float &t){

	vec3 p_to_c = ray.origin - center;
	float a = dot(ray.direction, ray.direction);
	float b = 2.0 * dot(p_to_c, ray.direction);
	float c = dot(p_to_c, p_to_c) - radius * radius;
	float discriminant = (b * b) - (4 * a * c); // b^2 - 4ac

	float t1 = (-b - sqrt(discriminant)) / (2*a);
	float t2 = (-b - sqrt(discriminant)) / (2*a);

	if(t1 < t2)
		t = t1;
	else
		t = t2;
	
	if(discriminant < 0 )
		return false;
	else
		return true;


}

vec3 Sphere :: get_normal(vec3 rayOrigin){

	return (rayOrigin - center) * (-1/(radius));
}