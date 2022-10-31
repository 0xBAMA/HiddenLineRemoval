# not-quite-a-demo-engine (aka NQADE /eŋkweɪd/)

# Setup
- Requires `libsdl2-dev` on Ubuntu.
- Make sure to recurse submodules to pull in FastNoise2, tracy, and hoard-of-bitfonts code - `git submodule update --init --recursive`, alternatively just run scripts/init.sh to do the same thing.
- Run scripts/build.sh in order to build on Ubuntu. It just automates the building of the executable and moves it to the root directory, as files are referenced relative to that location, and then deletes the build folder (if called with the "clean" option).



## Features
**Graphics Stuff**
- [SDL2](https://wiki.libsdl.org/) is used for windowing and input handling
- [GL3W](https://github.com/skaslev/gl3w) as an OpenGL loader
- OpenGL Debug callback for error/warning reporting
- [Dear ImGui Docking Branch](https://github.com/ocornut/imgui/tree/docking) + [ImGuiColorTextEdit](https://github.com/BalazsJako/ImGuiColorTextEdit)
	- Setup to spawn platform windows
- Fast text renderer overlay applied as a post process, currently used for reporting of main loop timing but I eventually want to extend it to a drop down in engine terminal a la [guake](http://guake-project.org/), renderer is loosely based on / inspired by [this implementation](https://jmickle66666666.github.io/blog/techart/2019/12/18/bitmap-font-renderer.html).
- WIP software rasterizer, SoftRast


**Noise**
- Diamond-Square algorithm heightmap generation
- [FastNoise2](https://github.com/Auburn/FastNoise2) - flexible, powerful, very fast noise generation
- 4 channel blue noise texture from [Christoph Peters](http://momentsingraphics.de/BlueNoise.html)


**Utilities**
- CMake build setup
- [GLM](http://glm.g-truc.net/0.9.8/api/index.html) for vector and matrix types
- Image wrapper for loading/saving/resizing/etc of images ( 4 channel, 8-bit default, WIP floating point support - fpng, stb, and lodepng backends )
- JSON parsing using [nlohmann's single header implementation](https://github.com/nlohmann/json)
- XML parsing using [TinyXML2](https://tinyxml2.docsforge.com/)
- TinyOBJLoader for loading of Wavefront .OBJ 3D model files
- [Brent Werness' Voxel Automata Terrain ( VAT )](https://bitbucket.org/BWerness/voxel-automata-terrain/src/master/), converted from processing to C++
	- BigInt library required by the VAT implementation ( to replace java/processing BigInt )
- [Tracy](https://github.com/wolfpld/tracy) profiler integration ( client mode )
