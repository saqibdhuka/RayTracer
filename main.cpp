// #include <GL/glew.h>
#include <GL/gl.h>
#include <GL/freeglut.h>
#include <GL/glut.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "Angel.h"
#include "ray.h"

using namespace std;

#define WINDOW 400
#define FOVY 45



void render(){
	for(int x=0; x < WINDOW; x++){
		for (int y = 0; y < WINDOW; y++)
		{
			vec3 orig(0.0, 0.0, 0.0);
			vec3 dir(x, y, -1.0);
			Ray ray(orig, dir);
		}
	}
}

void init(){
	
    glViewport(0, 0, WINDOW, WINDOW);
    glClearColor(0.0, 0.0, 0.0, 0.0);
	glEnable(GL_DEPTH_TEST);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective((float) FOVY, 1.0f, 0.1f, 100.0f);
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
	
	glTranslatef(0.0, 0.0, -5.0);
	glBegin(GL_TRIANGLES);
	glColor3f(1.0,0.0,0.0);
	    glVertex3f(0.0, 0.0, 0.0);
	    glVertex3f(1.0, 0.0, 0.0);
	    glVertex3f(0.5, 1.0, 0.0);
	glEnd();
	// glPopMatrix();

	glutSwapBuffers();

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

	glutMainLoop();
	return 0;
}
