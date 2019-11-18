#ifndef RAY_H
#define RAY_H

using namespace std;

extern *objectArray[];
class Ray
{

	public:

		Ray(vec3 o, vec3 d);
		Ray();

		vec3 origin;
		vec3 direction;
	
};

#endif