#include <GL/glew.h>
#include <GL/gl.h>
#include <GL/freeglut.h>
#include <GL/glut.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "Angel.h"
#include "ray.h"
#include "plane.h"
#include "sphere.h"
#include "light.h"


#define WINDOW 400
#define FOVY 45
#define zNear 0.1f
#define zFar 100.0f

float t = 100.0; //t for the use of Ray

const float LIMIT = tan(FOVY/2);

Sphere sph_objects[100];
Plane pla_objects[100];

Light light_object[100];

vec3 shading(Ray ray, Sphere sphere){

	vec3 color(1.0, 0.0, 1.0);

	if(sphere.insersect(ray, t)){
		vec3 V = ray.direction;
		vec3 P = ray.origin + V * t;
		vec3 N = sphere.get_normal(P);

		float ratio = dot(N,V);

		color = sphere.getkdColor() * (ratio *0.5);
	}

	return color;


}

void render(){
	for(float x=-LIMIT; x < LIMIT; x++){
		for (float y = -LIMIT; y < LIMIT; y++)
		{
			vec3 orig(0.0, 0.0, 0.0);
			vec3 dir(x, y, -1.0);
			Ray ray(orig, dir);
			ray.shootRay(t);
		}
	}
}

void init(){
	
    glViewport(0, 0, WINDOW, WINDOW);
    glClearColor(0.0, 0.0, 0.0, 0.0);
	glEnable(GL_DEPTH_TEST);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective((float) FOVY, 1.0f, zNear, zFar);
	glMatrixMode(GL_MODELVIEW);
	// gluLookAt(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
	// glOrtho(-100.0, 100.0, -100.0,100.0,-5.0,20.0);
	glLoadIdentity();
}


void callbackDisplay() {

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glEnable(GL_DEPTH);
	glMatrixMode(GL_MODELVIEW);

	// glPushMatrix();
		
	// glPopMatrix();

	

}



int main(int argc, char **argv) {

	// init GLUT and create Window
	glutInit(&argc, argv);
	glutInitDisplayMode( GLUT_SINGLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowPosition(1080, 450);
	glutInitWindowSize(WINDOW, WINDOW);
	glutCreateWindow("Ray Tracer");
	
	init();
	// register callbacks
	glutDisplayFunc(callbackDisplay);

	

	render();

	// glutSpecialFunc(specialKeyboard);
 	//glutKeyboardFunc(keyboard);
    // glutMouseFunc( mouse );

	// enter GLUT event processing cycle
	glutSwapBuffers();
	glutMainLoop();
	return 0;
}
