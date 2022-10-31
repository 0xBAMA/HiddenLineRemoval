#ifndef INCLUDES
#define INCLUDES

#include <stdio.h>

// stl includes
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <ctime>
#include <deque>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <vector>

// iostream stuff
using std::cerr;
using std::cin;
using std::cout;
using std::endl;
using std::flush;
using std::string;
constexpr char newline = '\n';

// pi definition - definitely sufficient precision
constexpr double pi = 3.14159265358979323846;

// MSAA count - effects OpenGL geometry evaluation
constexpr int MSAACount = 1;

// vector math library GLM
#define GLM_FORCE_SWIZZLE
#define GLM_SWIZZLE_XYZW
#include "../GLM/glm.hpp"                  //general vector types
#include "../GLM/gtc/matrix_transform.hpp" // for glm::ortho
#include "../GLM/gtc/type_ptr.hpp"         //to send matricies gpu-side
#include "../GLM/gtx/rotate_vector.hpp"
#include "../GLM/gtx/transform.hpp"

// not sure as to the utility of this
#define GLX_GLEXT_PROTOTYPES

// convenience defines for GLM
using glm::vec2;
using glm::vec3;
using glm::vec4;
using glm::ivec2;
using glm::ivec3;
using glm::ivec4;
using glm::mat3;
using glm::mat4;

// OpenGL Loader
#include "../ImGUI/gl3w.h"

// GUI library (dear ImGUI)
#include "../ImGUI/TextEditor.h"
#include "../ImGUI/imgui.h"
#include "../ImGUI/imgui_impl_sdl.h"
#include "../ImGUI/imgui_impl_opengl3.h"

// SDL includes - windowing, gl context, system info
#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

// image load/save/resize/access/manipulation wrapper
#include "../ImageHandling/Image.h"

// simple std::chrono wrapper
#include "Timer.h"

// keyboard input handling
// #include "inputHandler.h"

// tracy profiler annotation
#include "../tracy/public/tracy/Tracy.hpp"

// font rendering header
#include "../fonts/fontRenderer/renderer.h"

// wrapper for TinyOBJLoader
#include "../TinyOBJLoader/objLoader.h"

// software rasterizer reimplementation
#include "../SoftRast/SoftRast.h"

// shader compilation wrapper
#include "shaders/lib/shaderWrapper.h"

// coloring of CLI output
#include "../fonts/colors.h"

// diamond square heightmap generation
#include "../noise/diamondSquare/diamond_square.h"

// Brent Werness' Voxel Automata Terrain
#include "../noise/VAT/VAT.h"

// more general noise solution
#include "../noise/FastNoise2/include/FastNoise/FastNoise.h"

// bringing the old perlin impl back
#include "../noise/perlin.h"

// Niels Lohmann - JSON for Modern C++
#include "../JSON/json.hpp"
using json = nlohmann::json;

// tinyXML2 XML parser
#include "../tinyXML2/tinyxml2.h"
using XMLDocument = tinyxml2::XMLDocument;

// #define WIDTH 640
// #define HEIGHT 480

// #define WIDTH 1280
// #define HEIGHT 720

#define WIDTH 1920
#define HEIGHT 1080

struct colorGradeParameters {
	int tonemapMode = 6;
	float gamma = 1.1f;
	float colorTemp = 6500.0f;
};

// Function to get color temperature from shadertoy user BeRo
// from the author:
//   Color temperature (sRGB) stuff
//   Copyright (C) 2014 by Benjamin 'BeRo' Rosseaux
//   Because the german law knows no public domain in the usual sense,
//   this code is licensed under the CC0 license
//   http://creativecommons.org/publicdomain/zero/1.0/
// Valid from 1000 to 40000 K (and additionally 0 for pure full white)
inline vec3 GetColorForTemperature ( double temperature ) {
	// Values from:
	// http://blenderartists.org/forum/showthread.php?270332-OSL-Goodness&p=2268693&viewfull=1#post2268693
	mat3 m =
		( temperature <= 6500.0 )
			? mat3( vec3( 0.0, -2902.1955373783176, -8257.7997278925690 ),
					vec3( 0.0, 1669.5803561666639, 2575.2827530017594 ),
					vec3( 1.0, 1.3302673723350029, 1.8993753891711275 ) )
			: mat3( vec3( 1745.0425298314172, 1216.6168361476490, -8257.7997278925690 ),
					vec3( -2666.3474220535695, -2173.1012343082230, 2575.2827530017594 ),
					vec3( 0.55995389139931482, 0.70381203140554553, 1.8993753891711275 ) );

	return glm::mix( glm::clamp( vec3( m[ 0 ] / ( vec3( glm::clamp( temperature, 1000.0, 40000.0 ) ) +
		m[ 1 ] ) + m[ 2 ] ), vec3( 0.0 ), vec3( 1.0 ) ), vec3( 1.0 ), glm::smoothstep( 1000.0, 0.0, temperature ) );
}

#endif
