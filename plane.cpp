#include "plane.h"
#include <GL/gl.h>
#include <GL/freeglut.h>
#include <GL/glut.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "Angel.h"


Plane :: Plane(){
	A = 0.0;
	B = 0.0;
	C = 0.0;
	D = 0.0;

	matPlane.kdColor.x = 0.0;
	matPlane.kdColor.y = 0.0;
	matPlane.kdColor.z = 0.0;

	matPlane.ksColor.x = 0.0;
	matPlane.ksColor.y = 0.0;
	matPlane.ksColor.z = 0.0;

	matPlane.specular = 0.0;
	matPlane.reflection = 0.0;
	matPlane.refraction = 0.0;
	matPlane.index_refraction = 1.5;
}

Plane :: Plane (float a, float b, float c, float d, Material newMat){

	A = a;
	B = b;
	C = c;
	D = d;


	matPlane.kdColor.x = newMat.kdColor.x;
	matPlane.kdColor.y = newMat.kdColor.y;
	matPlane.kdColor.z = newMat.kdColor.z;

	matPlane.ksColor.x = newMat.ksColor.x;
	matPlane.ksColor.y = newMat.ksColor.y;
	matPlane.ksColor.z = newMat.ksColor.z;

	matPlane.specular = newMat.specular;
	matPlane.reflection = newMat.reflection;
	matPlane.refraction = newMat.refraction;
	matPlane.index_refraction = newMat.index_refraction;

}

vec3 Plane :: plane_normal(){
	vec3 N;
	N.x = A;
	N.y = B;
	N.z = C;

	return N;

}


float Plane :: intersectPlane(Ray &ray){

	vec3 normal = plane_normal();
	vec3 point = ray.origin;
	vec3 dir= ray.direction;
	float t;

	t = (dot(normal,ray.origin) + D) / dot(normal,ray.direction);
	if(t > 0){
		return t;
	}else{
		return -1;
	}

}
