50.017 Final Project - Fluidsim
===============================

Final project for 50.017 Graphics and Visualization course in Singapore University of Technology and Design (SUTD).

[![50.017 Fluidsim](https://j.gifs.com/yxPRge.gif)](http://www.youtube.com/watch?v=gq1194HM_N0 "50.017 Fluidsim Video:")

<b>Features:</b>  
- Smooth Particle Hydrodynamics (SPH) OpenCL / C++
- Marching Cubes with isovalue and edge caching (~8x faster than default Paul Bourke's version)
- SSE2 operations
- GLSL based reflection and refraction

Requirements
------------

OSX and XCode.

Running
-------

Change directory into the `executable` folder.

    cd executable

You may need to set the permissions on the `Fluidsim` application to be executable.

    sudo chmod +x Fluidsim

Then run

    ./Fluidsim

Compiling
---------

Open "GnV FluidSim.xcodeproj" in xCode.

By default, the SPH portion of the simulation uses the CPU. 

To toggle between OpenCL and CPU mode, go to fluidSystem.h, and comment/uncomment the following:

    #define OPENCL_SPH
