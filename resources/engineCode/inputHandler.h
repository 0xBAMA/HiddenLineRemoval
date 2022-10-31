#ifndef INPUT_H
#define INPUT_H
#include "includes.h"

//this is the other way to handle the input - check the state of the keys every frame
	// key element here is that you can handle more than one at once - e.g. I can have
	// simultaneous input from left and up, where the SDL_PollEvent method did not make
	// these things available in the same scope - it creates a massive spew of input
	// events when compared to the SDL_PollEvent method, so sometimes, that'll need to
	// be considered ( framerate dependency, cooldown kind of logic )
//==============================================================================
// list of input indices/names: https://wiki.libsdl.org/SDL_Scancode
//==============================================================================

// this goes in the event loop - was kind of a bust, not working correctly yet
// keyboard.Update();
// if ( keyboard.GetKeystate( SDL_SCANCODE_RIGHT ) == keyState::KEYDOWN )
// 	cout << "Right Key Pressed" << endl << flush;
// if ( keyboard.GetKeystate( SDL_SCANCODE_LEFT ) == keyState::KEYDOWN )
// 	cout << "Left Key Pressed" << endl << flush;
// if ( keyboard.GetKeystate( SDL_SCANCODE_UP ) == keyState::KEYDOWN )
// 	cout << "Up Key Pressed" << endl << flush;
// if ( keyboard.GetKeystate( SDL_SCANCODE_DOWN ) == keyState::KEYDOWN )
// 	cout << "Down Key Pressed" << endl << flush;

enum class keyState {
	KEYDOWN,
	HELD,
	KEYUP,
	UNPRESSED
};

// this is not really working as of yet, tbd
class inputHandler {
public:
	inputHandler () {
		const uint8_t * state = SDL_GetKeyboardState( &numKeys ); (void)state;
		previousState = ( uint8_t * ) malloc( numKeys );
		currentState = ( uint8_t * ) malloc( numKeys );
		Update(); // some ineffiency here, deleting previousState immediately
	} // but it's only once on startup, so probably not a huge concern

	void Update () {
		// currentState becomes previousState
		delete( previousState );
		previousState = currentState;
		// get new state and memcpy to new currentState array
		SDL_PumpEvents();
		const uint8_t * state = SDL_GetKeyboardState( NULL );
		currentState = ( uint8_t * ) malloc( numKeys );
		memcpy( ( void * ) state, currentState, numKeys );
	}

	keyState GetKeystate ( SDL_Scancode key ) {
		if ( currentState[ key ] && !previousState[ key ] ) {
			return keyState::KEYDOWN;
		} else if ( currentState[ key ] && previousState[ key ] ) {
			return keyState::HELD;
		} else if ( !currentState[ key ] && previousState[ key ] ) {
			return keyState::KEYUP;
		} else {
			return keyState::UNPRESSED;
		}
	}

	// TODO: mouse button handling with SDL_GetMouseState()

	int numKeys;
private:
	uint8_t * previousState	= nullptr;
	uint8_t * currentState	= nullptr;
};

#endif
