#include <GL/gl.h>
#include <GL/freeglut.h>
#include <GL/glut.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "Angel.h"



struct Material{
		vec3 kdColor;
		vec3 ksColor;
		int specular;//q
		float reflection;//kr
		float refraction;//kt
		float index_refraction;

		Material(){}

		Material(vec3 kd, vec3 ks, int q, float refl, float refr, float index) :
		kdColor(kd), ksColor(ks), specular(q), reflection(refl), refraction(refr), index_refraction(index){}

};
