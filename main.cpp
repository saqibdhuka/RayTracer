#include <GL/gl.h>
#include <GL/freeglut.h>
#include <GL/glut.h>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <time.h>
#include <iostream>
#include <array>
#include "Angel.h"
#include "ray.h"
#include "plane.h"
#include "sphere.h"
#include "light.h"

#ifndef MATERIALSTRUCTURE_H
#define MATERIALSTRUCTURE_H

#include "materialStructure.h"

#endif



#define WINDOW 400
#define FOVY 45
#define zNear 0.1f
#define zFar 100.0f




float tS; //t for the use of Ray
float t;
float tP; //t for the use of Ray

float xN, yN;
int rayDepth=4;
vec3 backColor(0.0, 0.0, 0.0);

const float LIMIT = tan(DegreesToRadians *FOVY/2)/WINDOW;

int countS =0;
int countP =0;
float tArray[3] = {500, 500, 500};//Array of t

Sphere sph_objects[3];
Plane pla_objects[1];

Light light_object[4];


void setRayDepth(int depth){
	rayDepth = depth;
}

vec3 vecPow(vec3 v, float power){

	vec3 ans(v.x, v.y, v.z);

	for(int i =2; i < power; i++){
		ans *= v;
	}

	return ans;
}

vec3 Ir(0,0,0);
vec3 It(0,0,0);
vec3 I_plane = backColor; //Light Intesity to calculate

vec3 tracePlane(Ray ray, Plane plane, Light light, int depth){ //Work on this

	if(depth == 0){
		return (I_plane);
		// return normalize(I);
	}else{
		if(yN > yN/2){
			return backColor;
		}

		vec3 kd = plane.matPlane.kdColor;
		vec3 Ii = light.colorIntensity;
		// float ambientStrength = 0.75f;
		vec3 ks = plane.matPlane.ksColor;
		// vec3 Ia =  Ii *ks;
		float  kr = plane.matPlane.reflection;
		vec3 Ia =  Ii *kr;
		float kt = plane.matPlane.refraction;
		int q = plane.matPlane.specular;
		float S =1;
		vec3 intersectionPoint = (ray.origin + (t * ray.direction));
		vec3 V = normalize(ray.origin - intersectionPoint);
		vec3 N = normalize(plane.plane_normal());
		vec3 Li = normalize(intersectionPoint - light_object[0].origin);
		vec3 R = normalize(((2*N)*(dot(N,V))) - V);
		float rv = dot(R,V);
		if(rv < 0){
			rv =0;
		}
		float tmp = pow(rv, q);
		vec3 ik = Ia * kd;
		vec3 ktemp = ks * tmp;
		float knl = dot((kd*N),Li);
		if(knl < 0){
			knl = 0;
		}
		vec3 sum = S*Ii*(vec3(knl + ktemp.x, knl + ktemp.y, knl + ktemp.z));
		// vec3 Ir = 0;
		// vec3 It = 0;


		// I_planes = ik + (sum);

		I_plane = ik + (Ir * kr) + (It * kt) + (sum);
		vec3 Ir = normalize(I_plane);//Reflected Ray Direction
		vec3 It = normalize(I_plane); //Refracted Ray Direction
		vec3 incidentRayVector = intersectionPoint - ray.origin;
		// vec3 reflectDir = ((2*N)*(dot(N,V))) - V;
		float tmpDiv = length(N) / length(incidentRayVector);
		float theta1 = acos(DegreesToRadians * tmpDiv);
		float refraction_index_air = 1.00029;
		float ratio = refraction_index_air / plane.matPlane.index_refraction; //n1/n2
		float c1 = cos(DegreesToRadians*theta1);
		float c2 = sqrt(1-(pow(ratio,2)*(pow(sin(DegreesToRadians*theta1), 2))));
		// float temp = ratio * sin(DegreesToRadians * theta1);
		// float theta2 = asin(temp);
		vec3 refractDir = ratio * incidentRayVector + (ratio*c1 - c2)*N;
		Ray ray2(intersectionPoint, R);
		Ray ray3(intersectionPoint, refractDir);

		return tracePlane(ray2, plane, light, depth-1);
		return tracePlane(ray3, plane, light, depth-1);
		// for(int i = depth; i >= 1; i--){
		// 	I_plane = ik + (Ir * kr) + (It * kt) + (sum);
		// 	vec3 Ir = I_plane;//Reflected Ray Direction
		// 	vec3 It = I_plane; //Refracted Ray Direction
		// }
		// I_plane+=0.3;
		// return normalize(I_plane);
		// return tracePlane(ray, plane, light, depth-1);

	}

}

vec3 I = backColor; //Light Intesity to calculate


vec3 trace(Ray ray, Sphere sphere, Light light, int depth){ //Sphere Trace


	if(depth == 0){
		return (I);
		// printf("Ix: %f,  Iy: %f,  Iz: %f\n", I.x, I.y, I.z);
		// return normalize(I);
	}else{
		vec3 kd = sphere.mat.kdColor;
		vec3 Ii = light.colorIntensity;
		// float ambientStrength = 0.75f;
		vec3 ks = sphere.mat.ksColor;
		// vec3 Ia =  Ii *ks;
		float  kr = sphere.mat.reflection;
		float kt = sphere.mat.refraction;
		vec3 Ia =  Ii *kr;
		int q = sphere.mat.specular;
		float S =1;
		vec3 intersectionPoint = (ray.origin + (t * ray.direction));
		vec3 V = normalize(ray.origin - intersectionPoint);
		vec3 N = normalize(sphere.get_normal(intersectionPoint));
		vec3 Li = normalize(intersectionPoint - light_object[0].origin);
		vec3 R = normalize(((2*N)*(dot(N,V))) - V);
		float rv = dot(R,V);
		if(rv < 0){
			rv =0;
		}
		float tmp = pow(rv, q);
		vec3 ik = Ia * kd;
		vec3 ktemp = ks * tmp;
		float knl = dot(kd*N,Li);
		if(knl < 0){
			knl = 0;
		}
		vec3 sum = S*Ii*(vec3(knl + ktemp.x, knl + ktemp.y, knl + ktemp.z));
		// vec3 Ir = 0;
		// vec3 It = 0;


		// I_planes = ik + (sum);

		I = ik + (Ir * kr) + (It * kt) + (sum);
		vec3 Ir = normalize(I);//Reflected Ray Direction
		vec3 It = normalize(I); //Refracted Ray Direction
		vec3 incidentRayVector = intersectionPoint - ray.origin;
		// vec3 reflectDir = ((2*N)*(dot(N,V))) - V;
		float tmpDiv = length(N) / length(incidentRayVector);
		float theta1 = acos(DegreesToRadians * tmpDiv);
		float refraction_index_air = 1.00029;
		float ratio = refraction_index_air / sphere.mat.index_refraction; //n1/n2
		float c1 = cos(DegreesToRadians*theta1);
		float c2 = sqrt(1-(pow(ratio,2)*(pow(sin(DegreesToRadians*theta1), 2))));
		// float temp = ratio * sin(DegreesToRadians * theta1);
		// float theta2 = asin(temp);
		vec3 refractDir = ratio * incidentRayVector + (ratio*c1 - c2)*N;
		Ray ray2(intersectionPoint, R);
		Ray ray3(intersectionPoint, refractDir);

		return trace(ray2, sphere, light, depth-1);
		return trace(ray3, sphere, light, depth-1);
		// for(int i = depth; i >= 1; i--){
		// 	I = ik + (Ir * kr) + (It * kt) + (sum);
		// 	vec3 Ir = I;//Reflected Ray Direction
		// 	vec3 It = I; //Refracted Ray Direction
		// }


		// printf("Ix: %f,  Iy: %f,  Iz: %f\n", I.x, I.y, I.z);


		// return normalize(I);
		// return trace(ray, sphere, light, depth-1);
	}

}

vec2 minT(Ray ray){

	vec2 min(10000,-1);//min[MinTValue, Index]
	for(int i=0; i < 3 ;i++){
		tArray[i] = sph_objects[i].intersect(ray);
		if(tArray[i] != -1 && tArray < min){
			min.x = tArray[i];
			min.y = i;
		}
		// float tPlane = pla_objects[0].intersectPlane(ray);
		// if(tPlane < tArray[i] && tPlane != -1){
		// 	min.x = tPlane;
		// 	min.y = -1;
		// }

	}
	return min;

}

void render(){

	float scale = 2;

	for(int x=-WINDOW/2; x < WINDOW/2; x++){
		 xN=x * LIMIT*scale;
		// printf("xN: %f\n", xN);
		for (int y = -WINDOW/2; y < WINDOW/2; y++)
		{
			 yN=y * LIMIT*scale;			// glBegin(GL_POINTS);
			// 	glVertex3f(xN, yN, -1.0);
			// glEnd();

			vec3 orig(0.0, 0.0, 0.0);
			vec3 dir(xN, yN, -1.0);
			Ray ray(orig, dir);

			vec2 T = minT(ray);
			t = T.x;
			// if(t!=-1)
				// printf("t: %f\n",t);
			int index = T.y;
			// printf("Index: %d\n",index);
			// t = sph_objects[0].intersect(ray);
			if(index != -1){
					countS++;
					vec3 colorLight = trace(ray, sph_objects[index], light_object[0], rayDepth);
					glColor3f(colorLight.x, colorLight.y, colorLight.z);
					// printf("  R: %f,  G: %f,  B: %f\n", colorLight.x, colorLight.y, colorLight.z);
					glBegin(GL_POINTS);
						glVertex3f(ray.origin.x + t*ray.direction.x,
							ray.origin.y + t*ray.direction.y,
							ray.origin.z + t*ray.direction.z);
					glEnd();
					// glutSwapBuffers();

			}
			else if(index <0){
				t = pla_objects[0].intersectPlane(ray);

				countP++;
				vec3 colorLight = tracePlane(ray, pla_objects[0], light_object[0], rayDepth);
				glColor3f(colorLight.x, colorLight.y, colorLight.z);

				glBegin(GL_POINTS);
					glVertex3f(ray.origin.x + t*ray.direction.x,
						ray.origin.y + t*ray.direction.y ,
						ray.origin.z + t*ray.direction.z );
				glEnd();

			}
		}
	}
	printf("countS: %d\n", countS);
	printf("countP: %d\n", countP);
}

void init(){

  glViewport(0, 0, WINDOW, WINDOW);

  glClearColor(backColor.x, backColor.y, backColor.z, 0.0);
	glEnable(GL_DEPTH_TEST);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective((float) FOVY, 1.0f, zNear, zFar);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 1.0, 0.0);
}


void callbackDisplay() {

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glEnable(GL_DEPTH);
	glMatrixMode(GL_MODELVIEW);



	vec3 center_ball1(0.0, -0.2, -10.0);
	vec3 center_ball2(1.0,0.0,-12.0);
	vec3 center_ball3(3.0, 0.3, -15.0);

	float radius_ball1 = 0.75;
	float radius_ball2 = 1.0;
	float radius_ball3 = 1.5;

	//Setting up the material
	Material Yellow = Material(vec3(0.93, 1.0, 0.0), vec3(1.0, 1.0, 1.0),128, 1.0f, 0.0f, 1.52);

	Material DarkRed = Material(vec3(0.502, 0.0353, 0.0353), vec3(1.0, 1.0, 1.0),80, 1.0f, 0.0f, 1.33);

	Material Purple = Material(vec3(0.4314, 0.0589, 0.6196), vec3(1.0, 1.0, 1.0),50, 1.0f, 0.0f, 1.5);

	Sphere ball3(center_ball3, radius_ball3, Purple);
	Sphere ball2(center_ball2, radius_ball2, DarkRed);
	Sphere ball1(center_ball1, radius_ball1, Yellow);

	sph_objects[0] = ball1;
	sph_objects[1] = ball2;
	sph_objects[2] = ball3;

	//Setting up the Plane
	Material PlaneMaterial = Material(vec3(0.0, 1.0, 1.0), vec3(1.0, 1.0, 1.0), 5.0, 1.0f, 0.0f, 1.31);
	Plane plane(0, 0, 1, -5, PlaneMaterial);
	pla_objects[0] = plane;


	//Setting up the Light
	vec3 lightPos(0.0, 0.0, -1.0);
	vec3 lightColor(0.5, 0.5, 0.5);//Ambient Light Color
	Light light(lightPos, lightColor);
	light_object[0] = light;

	render();
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
	setRayDepth(4);
	glutDisplayFunc(callbackDisplay);

	// enter GLUT event processing cycle
	// glutSwapBuffers();
	glutMainLoop();
	return 0;
}
