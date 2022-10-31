#ifndef OBJLOAD
#define OBJLOAD

#include "tiny_obj_loader.h"

class objLoader {
public:
	objLoader(){}
	objLoader( std::string fileName ){
		LoadOBJ( fileName );
	}

	void LoadOBJ( std::string fileName );

	// OBJ data ( per mesh )
	// this may vary in length
	std::vector< vec4 > vertices;
	std::vector< vec3 > normals;
	std::vector< vec3 > texcoords;

	// these should all be the same length, the number of triangles
	std::vector< ivec3 > triangleIndices;
	std::vector< ivec3 > normalIndices;
	std::vector< ivec3 > texcoordIndices;
};

#endif
