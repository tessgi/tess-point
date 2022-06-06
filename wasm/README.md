# Webassembly version of tess-point
This is based off the original c code version that has been compiled with
emscripten to make available to web browsers. This version along with
the custom tailored javascript binding is available on my
<a href="https://eta-earth.org/tess_play.html">personal webpage</a>.

The build uses emscripten and cmake. In order to build one needs to install
the emscripten sdk along with its prerequisites. In my personal experience using
a mac I needed Xcode, git, and cmake. On the M1 processor I ran into some troubles
with node.js version that came along with the emscripten sdk. I had to
install the arm64 node.js version first then install emscripten, so it used my local
version of node.js rather than bring another version. As another gotcha
for M1 laptops, the suggested Anaconda python for mac that comes up by default
is not the M1 (arm64) version that one needs. One has to select other versions
for mac os to find the anaconda arm64 version.

Prebuilt javascript and wasm are available in the public directory. If you want
to rebuild them, then follow a typical cmake build process.
The cmake instructions are in CMakeLists.txt, and emscripten takes care
of 'properly' calling emcc as the build target. From inside the top wasm directory,
issue the following commands to build.

`emcmake cmake -S . -B build`

`cd build`

`make`

I then move the resulting tess_stars2px.js and tess_stars2px.wasm into the
public directory by issuing `make install`

### TODOS:
1. Convert input and output to javascript arraybuffers
2. Support multi target file upload by figuring out how to load a local file containing coordinate lists with javascript
3. Allow output to be saved to textfile
