#include <GL/glew.h>
#include <GL/gl.h>
#include <GL/freeglut.h>
#include <GL/glut.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>




struct Material{
		vec3 kdColor;
		vec3 ksColor;	
		int specular;//q
		float reflection;//kr
		float refraction;//kt
};