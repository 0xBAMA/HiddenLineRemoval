#ifndef SOFTRAST
#define SOFTRAST

#include "../TinyOBJLoader/objLoader.h"
#include "../engineCode/includes.h"

struct triangle {
	glm::vec3 p0, p1, p2;
	glm::vec3 tc0, tc1, tc2;
	glm::vec3 n0, n1, n2;
};

// helper functions
static const float RemapRange ( const float value, const float iMin, const float iMax, const float oMin, const float oMax ) {
	return ( oMin + ( ( oMax - oMin ) / ( iMax - iMin ) ) * ( value - iMin ) );
}

static const rgba RGBAFromVec4( vec4 color ) {
	rgba temp;
	temp.r = uint8_t( RemapRange( color.r, 0.0f, 1.0f, 0.0f, 255.0f ) );
	temp.g = uint8_t( RemapRange( color.g, 0.0f, 1.0f, 0.0f, 255.0f ) );
	temp.b = uint8_t( RemapRange( color.b, 0.0f, 1.0f, 0.0f, 255.0f ) );
	temp.a = uint8_t( RemapRange( color.a, 0.0f, 1.0f, 0.0f, 255.0f ) );
	return temp;
}

static const vec3 BarycentricCoords ( vec3 p0, vec3 p1, vec3 p2, vec3 P ) {
	vec3 s[ 2 ];
	for ( int i = 2; i--; ) {
		s[ i ][ 0 ] = p2[ i ] - p0[ i ];
		s[ i ][ 1 ] = p1[ i ] - p0[ i ];
		s[ i ][ 2 ] = p0[ i ] - P[ i ];
	}
	vec3 u = glm::cross( s[ 0 ], s[ 1 ] );
	if ( std::abs( u[ 2 ] ) > 1e-2f )	 	// If u[ 2 ] is zero then triangle ABC is degenerate
		return vec3( 1.0f - ( u.x + u.y ) / u.z, u.y / u.z, u.x / u.z );
	return vec3( -1.0f, 1.0f, 1.0f );	// in this case generate negative coordinates, it will be thrown away by the rasterizer
}

static const mat3 rotation( vec3 a, float angle ) { //thanks to Neil Mendoza via http://www.neilmendoza.com/glsl-rotation-about-an-arbitrary-axis/
	a = glm::normalize( a ); //a is the axis
	float s = sin( angle );
	float c = cos( angle );
	float oc = 1.0f - c;
	return mat3(oc * a.x * a.x + c,         oc * a.x * a.y - a.z * s,  oc * a.z * a.x + a.y * s,
				oc * a.x * a.y + a.z * s,   oc * a.y * a.y + c,        oc * a.y * a.z - a.x * s,
				oc * a.z * a.x - a.y * s,   oc * a.y * a.z + a.x * s,  oc * a.z * a.z + c );
}

// Plans:
	// something to wrap texture reference, with or without interpolation - start with no interp for now
	// DrawModel, using TinyOBJLoader wrapper + transform
	// wrapper for writing a pixel's color and depth values, and optionally alpha blending
	// helper function for depth testing? it's already pretty short... tbd
	// object for holding triangle parameters, to simplify passing - need positions, need texcoords, need normals

class SoftRast {
public:
	SoftRast( uint32_t x = 0, uint32_t y = 0 ) : width( x ), height( y ) {
		Color = Image( x, y );
		Depth = ImageF( x, y );
		BlueNoise = Image( "resources/noise/blueNoise.png" ); // for sample jitter, write helper function to return some samples
		// init std::random generator as member variable, for picking blue noise sample point - then sweep along x or y to get low discrepancy sequence
	}

	vec4 BlueNoiseRef ( ivec2 loc ) {
		rgba value = BlueNoise.GetAtXY( loc.x % BlueNoise.width, loc.y % BlueNoise.height );
		return vec4( value.r / 255.0f, value.g / 255.0f, value.b / 255.0f, value.a / 255.0f ) - vec4( 0.5f );
	}

	Image currentTex;
	void LoadTex ( string texPath ) { currentTex = Image( texPath ); }
	vec4 TexRef ( vec2 texCoord ) {
		uint32_t x = uint32_t( texCoord.x * float( currentTex.width ) );
		uint32_t y = uint32_t( texCoord.y * float( currentTex.height ) );
		rgba val = currentTex.GetAtXY( x, y );
		vec4 returnVal = vec4( float( val.r ) / 255.0f, float( val.g ) / 255.0f, float( val.b ) / 255.0f, float( val.a ) / 255.0f );
		// cout << returnVal.x << " " << returnVal.y << " " << returnVal.z << " " << returnVal.w << newline;
		return returnVal;
	}


	const vec3 NDCToPixelCoords ( vec3 NDCCoord ) {
		return vec3(
			RemapRange( NDCCoord.x, -1.0f, 1.0f, 0.0f, float( width - 1.0f ) ),
			RemapRange( NDCCoord.y, -1.0f, 1.0f, 0.0f, float( height - 1.0f ) ),
			NDCCoord.z
		);
	}

	// draw dot - draw smooth? tbd, some kind of gaussian distribution, how do you handle depth?
	void DrawDot ( vec3 position, vec4 color ) {
		position = NDCToPixelCoords( position );
		// TODO: support for alpha blending, based on existing buffer color + input color
		vec2 positionXY = vec2( position.x, position.y );
		if ( glm::clamp( positionXY, vec2( 0.0f ), vec2( width, height ) ) == positionXY && // point is on screen
			Depth.GetAtXY( uint32_t( position.x ), uint32_t( position.y ) ).r > position.z ) { // depth testing
			Color.SetAtXY( position.x, position.y, RGBAFromVec4( color ) );
			Depth.SetAtXY( position.x, position.y, { position.z, 0.0f, 0.0f, 0.0f } );
		}
	}

	// draw line
	void DrawLine ( vec3 p0, vec3 p1, vec4 color ) {
		int x0 = int( RemapRange( p0.x, -1.0f, 1.0f, 0.0f, float( width ) - 1.0f ) );
		int y0 = int( RemapRange( p0.y, -1.0f, 1.0f, 0.0f, float( height ) - 1.0f ) );
		int x1 = int( RemapRange( p1.x, -1.0f, 1.0f, 0.0f, float( width ) - 1.0f ) );
		int y1 = int( RemapRange( p1.y, -1.0f, 1.0f, 0.0f, float( height ) - 1.0f ) );
		float z0 = p0.z;
		float z1 = p1.z;
		bool steep = false;
		if ( std::abs( x0 - x1 ) < std::abs( y0 - y1 ) ) {
			std::swap( x0, y0 );
			std::swap( x1, y1 );
			steep = true;
		}
		if ( x0 > x1 ) {
			std::swap( x0, x1 );
			std::swap( y0, y1 );
			std::swap( z0, z1 );
		}
		int dx = x1 - x0;
		int dy = y1 - y0;
		int derror2 = std::abs( dy ) * 2;
		int error2 = 0;
		int y = y0;
		for ( int x = x0; x <= x1; x++ ) {
			// interpolated depth value
			float depth = RemapRange( float( x ), float( x0 ), float( x1 ), z0, z1 );
			if ( steep ) {
				if ( Depth.GetAtXY( y, x ).r >= depth ) {
					Color.SetAtXY( y, x, RGBAFromVec4( color ) );
					Depth.SetAtXY( y, x, { depth, 0.0f, 0.0f, 0.0f } );
				}
			} else {
				if ( Depth.GetAtXY( x, y ).r >= depth ) {
					Color.SetAtXY( x, y, RGBAFromVec4( color ) );
					Depth.SetAtXY( x, y, { depth, 0.0f, 0.0f, 0.0f } );
				}
			}
			error2 += derror2;
			if ( error2 > dx ) {
				y += ( y1 > y0 ? 1 : -1 );
				error2 -= dx * 2;
			}
		}
	}

	// draw triangle
	// void DrawTriangle ( vec3 p0, vec3 p1, vec3 p2, vec4 color ) { // eventually extend to include texcoords + normals
	void DrawTriangle ( triangle t, vec4 color ) { // eventually extend to include texcoords + normals
		vec2 bboxmin(  std::numeric_limits< float >::max(),  std::numeric_limits< float >::max() );
		vec2 bboxmax( -std::numeric_limits< float >::max(), -std::numeric_limits< float >::max() );
		vec2 clamp( width - 1, height - 1 );

		// translate x, y of points into screen space
		t.p0 = NDCToPixelCoords( t.p0 );
		t.p1 = NDCToPixelCoords( t.p1 );
		t.p2 = NDCToPixelCoords( t.p2 );

		for ( int j = 0; j < 2; j++ ) {
			bboxmin[ j ] = std::max( 0.0f, std::min( bboxmin[ j ], t.p0[ j ] ) );
			bboxmax[ j ] = std::min( clamp[ j ], std::max( bboxmax[ j ], t.p0[ j ] ) );

			bboxmin[ j ] = std::max( 0.0f, std::min( bboxmin[ j ], t.p1[ j ] ) );
			bboxmax[ j ] = std::min( clamp[ j ], std::max( bboxmax[ j ], t.p1[ j ] ) );

			bboxmin[ j ] = std::max( 0.0f, std::min( bboxmin[ j ], t.p2[ j ] ) );
			bboxmax[ j ] = std::min( clamp[ j ], std::max( bboxmax[ j ], t.p2[ j ] ) );
		}

		const bool allowPrimitiveJitter = false;
		ivec2 eval;
		for ( eval.x = bboxmin.x; eval.x <= bboxmax.x; eval.x++ ) {
			for ( eval.y = bboxmin.y; eval.y <= bboxmax.y; eval.y++ ) {

				// for( n ) jittered samples? tbd, will need to do something to get an alpha value from the n samples
				vec4 jitter = allowPrimitiveJitter ? BlueNoiseRef( eval ) : vec4( 0.0f );
				vec3 bc = BarycentricCoords( t.p0, t.p1, t.p2, vec3( float( eval.x ) + jitter.x, float( eval.y ) + jitter.y, 0.0f ) );

				if ( bc.x < 0 || bc.y < 0 || bc.z < 0 ) continue; // any barycentric coord being negative means degenerate triangle or sample point outside triangle

				// if ( // interesting experiment, reject samples with certain ranges of the barycentric coords
				// 	( std::fmod( bc.x, 0.5f ) > 0.1618 && std::fmod( bc.y, 0.5f ) > 0.1618 ) ||
				// 	( std::fmod( bc.x, 0.5f ) > 0.1618 && std::fmod( bc.z, 0.5f ) > 0.1618 ) ||
				// 	( std::fmod( bc.z, 0.5f ) > 0.1618 && std::fmod( bc.y, 0.5f ) > 0.1618 )
				// ) continue;

				float depth = 0.0f; // barycentric interpolation of depth
				depth += bc.x * t.p0.z;
				depth += bc.y * t.p1.z;
				depth += bc.z * t.p2.z;

				vec2 texCoord = vec2( 0.0f );
				texCoord += bc.x * vec2( t.tc0.x, t.tc0.y );
				texCoord += bc.y * vec2( t.tc1.x, t.tc1.y );
				texCoord += bc.z * vec2( t.tc2.x, t.tc2.y );

				if ( Depth.GetAtXY( eval.x, eval.y ).r > depth ) {
					// will need to do this same barycentric interpolation of texcoords, normals, etc with bc

					// compute the color to write, texturing, etc, etc

					// write color - start with bc as color, 1.0 alpha - eventually will need to blend with existing color buffer value
					// Color.SetAtXY( eval.x, eval.y, { uint8_t( bc.x * 255.0f ), uint8_t( bc.y * 255.0f ), uint8_t( bc.z * 255.0f ), 255 } );
					// Color.SetAtXY( eval.x, eval.y, { uint8_t( texCoord.x * 255.0f ), uint8_t( texCoord.y * 255.0f ), 0, 255 } );

					vec4 texRef = TexRef( vec2( texCoord.x, 1.0f - texCoord.y ) );
					if ( texRef.a == 0.0f ) continue; // reject zero alpha samples - still need to blend

					Color.SetAtXY( eval.x, eval.y, { uint8_t( texRef.x * 255 ), uint8_t( texRef.y * 255 ), uint8_t( texRef.z * 255 ), uint8_t( texRef.w * 255 ) } );
					Depth.SetAtXY( eval.x, eval.y, { depth, 0.0f, 0.0f, 0.0f } );
				}
			}
		}
	}

	void DrawModel ( string modelPath, string texturePath, mat3 transform, vec3 offset ) {
		// passing in transform means we can scale, rotate, etc, and keep the interface simple
		objLoader o( modelPath );

		LoadTex( texturePath );

		// cout << "image loaded " << currentTex.width << " by " << currentTex.height << newline;

		for ( unsigned int i = 0; i < o.triangleIndices.size() - 80 /* exclude wheel */; i++ ) {
			vec4 p0 = o.vertices[ int( o.triangleIndices[ i ].x ) ];
			vec4 p1 = o.vertices[ int( o.triangleIndices[ i ].y ) ];
			vec4 p2 = o.vertices[ int( o.triangleIndices[ i ].z ) ];
			vec3 p0x( p0.x, p0.y, p0.z );
			vec3 p1x( p1.x, p1.y, p1.z );
			vec3 p2x( p2.x, p2.y, p2.z );

			p0x = transform * p0x;
			p1x = transform * p1x;
			p2x = transform * p2x;

			p0x += offset;
			p1x += offset;
			p2x += offset;

			triangle t;
			t.p0 = p0x;
			t.p1 = p1x;
			t.p2 = p2x;
			t.tc0 = o.texcoords[ int( o.texcoordIndices[ i ].x ) ];
			t.tc1 = o.texcoords[ int( o.texcoordIndices[ i ].y ) ];
			t.tc2 = o.texcoords[ int( o.texcoordIndices[ i ].z ) ];

			// DrawTriangle ( p0x, p1x, p2x, vec4( 1.0f ) );
			DrawTriangle ( t, vec4( 1.0f ) );
		}

		// wireframe
		// for ( unsigned int i = 0; i < o.triangleIndices.size() - 80 /* exclude wheel */; i++ ) {
		// 	vec4 p0 = o.vertices[ int( o.triangleIndices[ i ].x ) ];
		// 	vec4 p1 = o.vertices[ int( o.triangleIndices[ i ].y ) ];
		// 	vec4 p2 = o.vertices[ int( o.triangleIndices[ i ].z ) ];
		// 	vec3 p0x( p0.x, p0.y, p0.z );
		// 	vec3 p1x( p1.x, p1.y, p1.z );
		// 	vec3 p2x( p2.x, p2.y, p2.z );

		// 	const vec3 offset = vec3( 0.0, 0.0, -0.01 );
		// 	p0x = ( transform * p0x ) + offset;
		// 	p1x = ( transform * p1x ) + offset;
		// 	p2x = ( transform * p2x ) + offset;

		// 	triangle t;
		// 	t.p0 = p0x;
		// 	t.p1 = p1x;
		// 	t.p2 = p2x;
		// 	t.tc0 = o.texcoords[ int( o.texcoordIndices[ i ].x ) ];
		// 	t.tc1 = o.texcoords[ int( o.texcoordIndices[ i ].y ) ];
		// 	t.tc2 = o.texcoords[ int( o.texcoordIndices[ i ].z ) ];
		//
		// 	DrawLine ( p0x, p1x, vec4( 1.0f ) );
		// 	DrawLine ( p0x, p2x, vec4( 1.0f ) );
		// 	DrawLine ( p2x, p1x, vec4( 1.0f ) );
		// }

	}

	// dimensions
	uint32_t width = 0;
	uint32_t height = 0;

	// buffers
	Image Color;
	ImageF Depth;
	Image BlueNoise;
};

#endif
