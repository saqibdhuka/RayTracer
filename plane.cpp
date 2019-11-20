#include "plane.h"
#include <GL/glew.h>
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
}

Plane :: Plane (float a, float b, float c, float d, Material newMat){

	A = a;
	B = b;
	C = c;
	D = d;


	newMat.kdColor.x = matPlane.kdColor.x;
	newMat.kdColor.y = matPlane.kdColor.y;
	newMat.kdColor.z = matPlane.kdColor.z;

	newMat.ksColor.x = matPlane.ksColor.x;
	newMat.ksColor.y = matPlane.ksColor.y;
	newMat.ksColor.z = matPlane.ksColor.z;

	newMat.specular = matPlane.specular;
	newMat.reflection = matPlane.reflection;
	newMat.refraction = matPlane.refraction;
}


void Plane::draw(){
	
}



using namespace std;