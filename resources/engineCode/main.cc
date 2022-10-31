#include "engine.h"

int main ( int argc, char *argv[] ) {
	engine e;
	while( !e.MainLoop() );
	return 0;
}
