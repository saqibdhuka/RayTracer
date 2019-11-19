#ifndef SPHERE_H
#define SPHERE_H




class Sphere
{
public:


	struct Material{
		vec3 kdColor;
		vec3 ksColor;	
		float specular;
		float reflection;
		float refraction;
	};
	Sphere();
	Sphere(vec3 point, float rad, Material newMat);
	void draw();
	vec3 centrer;
	float radius;
	Material mat;
	
};


#endif