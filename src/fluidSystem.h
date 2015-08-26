#ifndef FLUIDSYSTEM_H
#define FLUIDSYSTEM_H

#include <vector>

using namespace std;

#include "helpers.h"
#include "marchingCubes.h"

// define it to activate OpenCL
//#define OPENCL_SPH

#ifdef OPENCL_SPH
#include <OpenCL/opencl.h>
#include "opencl_helpers.h"
#endif

class FluidSystem
{
public:

	FluidSystem(float width, float height, float depth);
	~FluidSystem();
	void update();
	void setGravDirection(const float4 &gravDirection);
	void draw();
	bool getDrawParticleMode();
	void setDrawParticleMode(bool mode);

private:
	
#ifdef OPENCL_SPH
	SourceReader source;

	cl_device_id     device;
	cl_context       clContext;
	cl_command_queue queue;
	cl_program       clProgram;
	cl_kernel        updateParticleKernel;
	cl_kernel        updateDensityKernel;
	cl_kernel        updateForceKernel;
	
	cl_mem cl_gravDirection;
	cl_mem cl_gridStarts;
	cl_mem cl_gridEnds;
	cl_mem cl_particleDensities;
	cl_mem cl_particlePositions;
	cl_mem cl_particleVelocities;
	cl_mem cl_particleForces;
#endif
	
	float4v verts, norms;
	mini_vec<Face> faces;
	
	// Struct of arrays instead of arrays of structs
	int *gridStarts;
	int *prevGridStarts;
	int *gridEnds;
	int *prevGridEnds;
	
	int *gridMakeStarts;
	int *gridMakeEnds;
	int *gridMakeParticleNexts;
	
	float *particleDensities;
	float4 *particlePositions;
	float4 *gridMakeParticlePositions;
	float4 *particleVelocities;
	float4 *gridMakeParticleVelocities;
	float4 *particleForces;
	float4 *gravDirection;
	

	// Physical constants...
	float MASS;
	float K;
	float HALF_K;
	float MU;
	float REST_DENSITY;
	float REST_DENSITY_TIMES_2;
	float SIGMA;
	float POINT_DAMPING;
	float H;
	float HH;
	float TIME_STEP;
	float GRAV_CONST;
	float SCALE;
	float COLLISION_ELASTICITY;
	
	float KERNEL_BASE;
	float PRESSURE_GRADIDENT_BASE;
	float VISCOSITY_LAPLACIAN_BASE;
	float GRADIENT_OR_LAPLACIAN_BASE;
	
	// Resolution related constants...
	float MC_RES;
	int   MC_MARGIN;
	float MC_THRESHOLD;
	float MC_POLY_THRESHOLD;
	
	size_t PARTICLE_COUNT;
	int GRID_WIDTH;
	int GRID_HEIGHT;
	int GRID_WIDTH_HEIGHT;
	int GRID_DEPTH;
	int GRID_SIZE;
	
	float WIDTH;
	float HEIGHT;
	float DEPTH;
	
	bool sphereInited;
	int sphereId;
	bool centerDrawing;
	bool drawParticlesMode;
	int drawParticlesTimer;
	int DRAW_PARTICLES_TIME_LIMIT;
	
	float4 drawingOffset;
	
	void updateGrid(); // To be done in CPU, as non-trivially parallel
	
	void updateDensities(); // To be done in OpenCL
	
	void updateForces();  // To be done in OpenCL
	
	void updateParticles();  // To be done in OpenCL

	float mcDensityAtPoint(const float4 &point);
	
	float mcDensityAtPointByParticle(const float4 &point, int p);
	
	void generateIsoSurface(float4v &verts, float4v &norms, mini_vec<Face> &faces);
		
	void drawParticle(int p);
};

#endif
