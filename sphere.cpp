#include "sphere.h"


using namespace std;



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
	
}