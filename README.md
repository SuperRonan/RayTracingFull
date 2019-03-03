# RayTracingFull

Real time Ray Tracing and global illumination engine.


--------Features--------

CPU Real time ray tracer in resolutions arround 300x300, and the resolution can be changed on the fly. 
You can fly into the scene using ZQSD and arrows.
The rendering is accelerated with a BVH (see Geometry/bvh_scene).
There are different rendering modes: "realistic" phong, normals, uv, acceleration strucure heat map and depth.

The pre calculation of the shadows is currently desactivated due to some little artifacts, a longer pre-computation time 
and of course less precise shadows. It can be turned on by uncommenting the line to launch it in the main. 
But it offers faster performances (and soft shadows thanks to the interpolation at the edges).  

There is the beggining of a global illumination mode, but it is still very noisy, biased, and wrong.
The result of the global illumination should be very different to the "basic" ray tracing.
The "basic" ray tracing uses only point lights for ligthing and the global illumination only samples from surface lights.

The project contains different scenes. You can select the scene by uncomenting its init function at the begining of the main.
The number of triangles of the scene is printed in the console. 

There is also a skybox (file skybox.png). It doesn't look that bad as long as you don't watch too much up or down. 

--------Instalation procedure--------

This is a Visual Studio 2017 project.

In the dependencies_ima folder, execute setup.bat to setup the dependencies.

You can launch the VS 2017 project with RayCasting_git/project2017/RayCasting.sln.

The projects works with x64 configuration and should also work with win32 (not sure).
Be sure to launch the release version, and launch it with Ctrl + F5 to avoid the debugger and improve performances.  



--------Important files--------

main.cpp:
- Select the scene, and other settings

Scene.h:
- Rendering of the scene (compute, sendRay, phong (different versions))
- Shadow pre computing (pre_compute_shadows, compute_triangle_shadow)

bvh_scene.h (inherits Scene.h)
- Building of the bvh: pre_compute
- Overload of the intersection function

turbo_bvh_scene.h:
- Justs builds the bvh faster with multi-threading

TriangleMultipleTexture.hpp:
- Class to represent multiple textures on a triangle (tried to optimize storage)



--------Authorship and prodution--------

This is a practical work project by Ronan Cailleau.
The base of the project is by Fabrice Lamarche. 
The course is by Fabrice Lamarche and Remi Cozot. 

This is a practical work project where I tried a lot of things, so the code is a dirty at some points.  
