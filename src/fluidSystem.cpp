#include "fluidSystem.h"
#include <cstdio>

#include <iostream>

using namespace std;

#define SQR(x)	((x)*(x))
#define CUBE(x)	((x)*(x)*(x))
#define POW6(x)	(CUBE(x)*CUBE(x))
#define POW9(x)	(POW6(x)*CUBE(x))

#define GRID_INDEX(i,j,k) (GRID_WIDTH_HEIGHT * k + GRID_WIDTH * j + i)
#define KERNEL(kernelDiff) KERNEL_BASE * CUBE(kernelDiff)
#define GRADIENT(d, kernelDiff) (-GRADIENT_OR_LAPLACIAN_BASE * SQR(kernelDiff) * d)
#define LAPLACIAN(kernelDiff) (GRADIENT_OR_LAPLACIAN_BASE * kernelDiff * (7.0f * rSq - 3.0f * HH))
#define PRESSURE_GRADIENT(d, r, kernelDiff2) \
((r<0.001f) ? 0.f : (PRESSURE_GRADIDENT_BASE * SQR(kernelDiff2) * d / r))
#define VISCOSITY_LAPLACIAN(kernelDiff2) (VISCOSITY_LAPLACIAN_BASE * kernelDiff2)

#define POSITION_GRID_INDEX(position) \
GRID_INDEX((int)position[0],\
(int)position[1],\
(int)position[2]);\

#define FOREACH_I_J_K(i,j,k) \
for (int k=0; k<GRID_DEPTH; ++k)\
for (int j=0; j<GRID_HEIGHT; ++j)\
for (int i=0; i<GRID_WIDTH; ++i)\

#define FOREACH_PARTICLE(p) \
for (int p=0; p<PARTICLE_COUNT; ++p)

#ifdef OPENCL_SPH
#define SET_C_AND_CL_FLOAT(var, value) var = value; source.defineFloat(#var, value);
#define SET_C_AND_CL(var, value) var = value; source.define(#var, value);
#else
#define SET_C_AND_CL_FLOAT(var, value) var = value;
#define SET_C_AND_CL(var, value) var = value;
#endif

inline float4 planeNorm(const float4& p0, const float4& p1, const float4& p2) {
	return cross(p1-p0, p2-p0).normalized();
}

#define DO_3x(c) c;c;c;
inline void calculateNorms(float4v &verts, float4v &norms, mini_vec<Face> &faces) {
	norms.clear();
	norms.resize(verts.size());
	for (int i=0; i<norms.size(); ++i) norms[i] = float4(0.f);
	for (int fi=0; fi<faces.size(); ++fi) {
		Face &f = faces[fi];
		float4 n = planeNorm(verts[f[0]], verts[f[1]], verts[f[2]]);
		int i=0;
		DO_3x(norms[f[i]] += n; ++i;)
	}
	for (int i=0; i<norms.size(); ++i) norms[i].normalize();
}

inline void drawFaces(float4v &verts, float4v &norms, mini_vec<Face> &faces) {
	glBegin(GL_TRIANGLES);
	for (int fi=0; fi<faces.size(); ++fi) {
		Face &f = faces[fi];
		glNormal3fv(norms[f[0]]); glVertex3fv(verts[f[0]]);
		glNormal3fv(norms[f[1]]); glVertex3fv(verts[f[1]]);
		glNormal3fv(norms[f[2]]); glVertex3fv(verts[f[2]]);
	}
	glEnd();
}


FluidSystem::FluidSystem(float width, float height, float depth)
{
#ifdef OPENCL_SPH
	source = SourceReader("fluid.cl");
#endif
	
	SET_C_AND_CL_FLOAT(MASS, 1.f);
	SET_C_AND_CL_FLOAT(GRAV_CONST, 150.f);
	SET_C_AND_CL_FLOAT(SCALE, .2f);
	SET_C_AND_CL_FLOAT(K, 1000.f); // Gas constant
	SET_C_AND_CL_FLOAT(HALF_K, K*.5f);
	SET_C_AND_CL_FLOAT(MU, .1f); // Mu for Fluid
	SET_C_AND_CL_FLOAT(REST_DENSITY, 1.2f);
	SET_C_AND_CL_FLOAT(REST_DENSITY_TIMES_2, REST_DENSITY*2.f);
	SET_C_AND_CL_FLOAT(SIGMA, 1.8f); // 1.5f
	SET_C_AND_CL_FLOAT(POINT_DAMPING, 2.0f);
	SET_C_AND_CL_FLOAT(H, 1.5f); // Kernel width
	SET_C_AND_CL_FLOAT(TIME_STEP, 0.01f);
	SET_C_AND_CL_FLOAT(HH, H*H);
	SET_C_AND_CL(PARTICLE_COUNT, 3000);
	SET_C_AND_CL(GRID_WIDTH, int(width / H) + 1);
	SET_C_AND_CL(GRID_HEIGHT, int(height / H) + 1);
	SET_C_AND_CL(GRID_DEPTH, int(depth / H) + 1);
	SET_C_AND_CL(GRID_WIDTH_HEIGHT, GRID_WIDTH * GRID_HEIGHT);
	SET_C_AND_CL_FLOAT(WIDTH, width);
	SET_C_AND_CL_FLOAT(HEIGHT, height);
	SET_C_AND_CL_FLOAT(DEPTH, depth);
	SET_C_AND_CL_FLOAT(COLLISION_ELASTICITY, 0.8f); // 1.f
	SET_C_AND_CL_FLOAT(KERNEL_BASE, 315.f  / (64.f * M_PI * POW9(H)));
	SET_C_AND_CL_FLOAT(GRADIENT_OR_LAPLACIAN_BASE, -945.f / (32.f * M_PI * POW9(H)));
	SET_C_AND_CL_FLOAT(VISCOSITY_LAPLACIAN_BASE, 45.0f  / (M_PI * POW6(H)));
	SET_C_AND_CL_FLOAT(PRESSURE_GRADIDENT_BASE, -VISCOSITY_LAPLACIAN_BASE);
	SET_C_AND_CL(GRID_SIZE, GRID_WIDTH * GRID_HEIGHT * GRID_DEPTH);
	
#ifdef OPENCL_SPH
	source.defineFloat("H_RECIP", 1.f/H);
	source.defineFloat("MASS_RECIP", 1.f/MASS);
#endif
	
	MC_RES = 0.5f;
	MC_MARGIN = 3;
	MC_THRESHOLD = 2.5f * MASS / KERNEL_BASE;
	MC_POLY_THRESHOLD = MASS * .6 / KERNEL_BASE;
	
	sphereInited = false;
	drawParticlesMode = false;
	drawParticlesTimer = 0;
	centerDrawing = true;
	DRAW_PARTICLES_TIME_LIMIT = 60 * 60; // in frames, so it's an approx
	
	particleDensities  = new float[PARTICLE_COUNT];
	particlePositions  = new float4[PARTICLE_COUNT];
	particleVelocities = new float4[PARTICLE_COUNT];
	particleForces     = new float4[PARTICLE_COUNT];
	gravDirection      = new float4[1];
	
	gridMakeParticlePositions  = new float4[PARTICLE_COUNT];
	gridMakeParticleVelocities = new float4[PARTICLE_COUNT];
	gridStarts                 = new int[GRID_SIZE+2];
	gridEnds                   = new int[GRID_SIZE+2];
	gridMakeStarts             = new int[GRID_SIZE+2];
	gridMakeEnds               = new int[GRID_SIZE+2];
	prevGridStarts             = new int[GRID_SIZE+2];
	prevGridEnds	               = new int[GRID_SIZE+2];
	gridMakeParticleNexts      = new int[PARTICLE_COUNT];
	
	for (int p=0; p<PARTICLE_COUNT; ++p) gridMakeParticleNexts[p] = -1;
	
	for (int g=0; g<GRID_SIZE; ++g) {
		gridStarts[g] =
		gridEnds[g] =
		gridMakeStarts[g] =
		gridMakeEnds[g] =
		prevGridStarts[g] =
		prevGridEnds[g] = -1;
	}
	
	for (int g=GRID_SIZE; g<GRID_SIZE+2; ++g) {
		gridStarts[g] =
		gridEnds[g] =
		gridMakeStarts[g] =
		gridMakeEnds[g] =
		prevGridStarts[g] =
		prevGridEnds[g] = (int) PARTICLE_COUNT - 1;
	}
	
	int p = 0;
	bool loop = true;
	while (loop) {
		for (int j=0; loop && j<height; ++j) {
			for (int k=0; loop && k<depth; ++k) {
				for (int i=0; loop && i<width; ++i) {
					float4 pos =  float4(i,j,k);
					
					particlePositions[p]  = pos;
					particleVelocities[p] = float4(0.f);
					particleForces[p]     = float4(0.f);
					particleDensities[p]  = 0;
					
					p++;
					
					if (p == PARTICLE_COUNT) {
						loop = false;
					}
					
				}
			}
		}
	}
	
	updateGrid(); // Take note... must update grid here too.
	
	if (centerDrawing) drawingOffset = -0.5 * float4(WIDTH, HEIGHT, DEPTH);
	else drawingOffset = 0;
	
#ifdef OPENCL_SPH
	CL_GET_GPU_CONTEXT_AND_COMMAND_QUEUE(device, clContext, queue);
	
	CL_CREATE_BUFFER(cl_gridStarts, clContext, CL_MEM_READ_ONLY, sizeof(int) * GRID_SIZE);
	CL_CREATE_BUFFER(cl_gridEnds, clContext, CL_MEM_READ_ONLY, sizeof(int) * GRID_SIZE);
	CL_CREATE_BUFFER(cl_particleDensities, clContext, CL_MEM_READ_WRITE, sizeof(float) * PARTICLE_COUNT);
	CL_CREATE_BUFFER(cl_particlePositions, clContext, CL_MEM_READ_WRITE, sizeof(float4) * PARTICLE_COUNT);
	CL_CREATE_BUFFER(cl_particleVelocities, clContext, CL_MEM_READ_WRITE, sizeof(float4) * PARTICLE_COUNT);
	CL_CREATE_BUFFER(cl_particleForces, clContext, CL_MEM_READ_WRITE, sizeof(float4) * PARTICLE_COUNT);
	CL_CREATE_BUFFER(cl_gravDirection, clContext, CL_MEM_READ_ONLY, sizeof(float4));
	
	const char *kernelSource = source;
	CL_CREATE_AND_BUILD_PROGRAM(clProgram, clContext, kernelSource);
	CL_CREATE_KERNEL(updateDensityKernel, clProgram, "updateDensity");
	CL_CREATE_KERNEL(updateForceKernel, clProgram, "updateForce");
	CL_CREATE_KERNEL(updateParticleKernel, clProgram, "updateParticle");
	
	CL_SET_KERNEL_ARG(updateDensityKernel, 0, sizeof(cl_mem), cl_gridStarts);
	CL_SET_KERNEL_ARG(updateDensityKernel, 1, sizeof(cl_mem), cl_gridEnds);
	CL_SET_KERNEL_ARG(updateDensityKernel, 2, sizeof(cl_mem), cl_particlePositions);
	CL_SET_KERNEL_ARG(updateDensityKernel, 3, sizeof(cl_mem), cl_particleDensities);
	
	CL_SET_KERNEL_ARG(updateForceKernel, 0, sizeof(cl_mem), cl_gridStarts);
	CL_SET_KERNEL_ARG(updateForceKernel, 1, sizeof(cl_mem), cl_gridEnds);
	CL_SET_KERNEL_ARG(updateForceKernel, 2, sizeof(cl_mem), cl_particlePositions);
	CL_SET_KERNEL_ARG(updateForceKernel, 3, sizeof(cl_mem), cl_particleVelocities);
	CL_SET_KERNEL_ARG(updateForceKernel, 4, sizeof(cl_mem), cl_particleForces);
	CL_SET_KERNEL_ARG(updateForceKernel, 5, sizeof(cl_mem), cl_particleDensities);
	CL_SET_KERNEL_ARG(updateForceKernel, 6, sizeof(cl_mem), cl_gravDirection);
	
	CL_SET_KERNEL_ARG(updateParticleKernel, 0, sizeof(cl_mem), cl_particlePositions);
	CL_SET_KERNEL_ARG(updateParticleKernel, 1, sizeof(cl_mem), cl_particleVelocities);
	CL_SET_KERNEL_ARG(updateParticleKernel, 2, sizeof(cl_mem), cl_particleForces);
	CL_SET_KERNEL_ARG(updateParticleKernel, 3, sizeof(cl_mem), cl_particleDensities);
	
	CL_ENQUEUE_WRITE_BUFFER(queue, cl_gridStarts, sizeof(int) * GRID_SIZE, gridStarts);
	CL_ENQUEUE_WRITE_BUFFER(queue, cl_gridEnds, sizeof(int) * GRID_SIZE, gridEnds);
	CL_ENQUEUE_WRITE_BUFFER(queue, cl_particlePositions, sizeof(float4) * PARTICLE_COUNT, particlePositions);
	CL_ENQUEUE_WRITE_BUFFER(queue, cl_particleVelocities, sizeof(float4) * PARTICLE_COUNT, particleVelocities);
	CL_ENQUEUE_WRITE_BUFFER(queue, cl_gravDirection, sizeof(float4), gravDirection);
#endif
}

FluidSystem::~FluidSystem()
{
#ifdef OPENCL_SPH
	CL_RELEASE(cl_gridStarts, cl_gridEnds, cl_particleDensities,
			   cl_particlePositions, cl_particleVelocities, cl_particleForces,
			   updateParticleKernel, updateDensityKernel, clProgram, queue, clContext);
#endif
	delete [] gridStarts;
	delete [] prevGridStarts;
	delete [] gridEnds;
	delete [] prevGridEnds;
	
	delete [] gridMakeStarts;
	delete [] gridMakeEnds;
	delete [] gridMakeParticleNexts;
	
	delete [] particleDensities;
	delete [] particlePositions;
	delete [] gridMakeParticlePositions;
	delete [] particleVelocities;
	delete [] gridMakeParticleVelocities;
	delete [] particleForces;
	delete [] gravDirection;
	
}

inline void FluidSystem::updateGrid()
{
	for (int g=0; g<GRID_SIZE; ++g) {
#ifdef OPENCL_SPH
		prevGridStarts[g] = gridStarts[g];
		prevGridEnds[g] = gridEnds[g];
#endif
		gridMakeStarts[g] = gridMakeEnds[g] = -1;
	}
	
	FOREACH_PARTICLE(p) {
		
		float4 positionIndexes = particlePositions[p] / H;
		
		int g = POSITION_GRID_INDEX(positionIndexes);
		
		if (gridMakeStarts[g] < 0) {
			gridMakeStarts[g] = p;
			gridMakeEnds[g] = p;
		} else {
			gridMakeParticleNexts[gridMakeEnds[g]] = p;
			gridMakeEnds[g] = p;
		}
		gridMakeParticleNexts[p] = -1;
	}
	
	int c = 0;
	for (int g=0; g<GRID_SIZE; ++g) {
		gridStarts[g] = c;
		for (int p=gridMakeStarts[g]; p>-1; p=gridMakeParticleNexts[p]) {
			gridMakeParticlePositions[c]  = particlePositions[p];
			gridMakeParticleVelocities[c] = particleVelocities[p];
			++c;
		}
		gridEnds[g] = c;
	}
	FOREACH_PARTICLE(p) {
		particleVelocities[p] = gridMakeParticleVelocities[p];
#ifdef OPENCL_SPH
		swap(particlePositions[p], gridMakeParticlePositions[p]);
#else
		particlePositions[p] = gridMakeParticlePositions[p];
#endif
	}
}

void FluidSystem::setGravDirection(const float4 &gravDirection)
{
	this->gravDirection[0] = gravDirection;
}

inline void FluidSystem::drawParticle(int particleId)
{
	float4 p = SCALE * (particlePositions[particleId] + drawingOffset);
	glTranslatef(+p[0], +p[1], +p[2]);
	glCallList(sphereId);
	glTranslatef(-p[0], -p[1], -p[2]);
	
}

#define FOREACH_P(p) \
int gw = GRID_WIDTH;\
int gwh = GRID_WIDTH_HEIGHT;\
int gwmm = GRID_WIDTH-1, ghmm = GRID_HEIGHT-1, gdmm = GRID_DEPTH-1;\
int p, pEnd, i, j, k, gridK, gridKEnd, gridJK, gridJKEnd, gridIJK, gridIJKEnd;\
for (k=0, gridK = 0, gridKEnd = GRID_DEPTH * gwh; gridK < gridKEnd; ++k, gridK += gwh)\
for (j=0, gridJK = gridK, gridJKEnd = gridK + gwh; gridJK < gridJKEnd; ++j, gridJK += gw)\
for (i=0, gridIJK = gridJK, gridIJKEnd = gridJK + gw; gridIJK < gridIJKEnd; ++i, ++gridIJK)\
for (p = gridStarts[gridIJK], pEnd = gridEnds[gridIJK]; p < pEnd; ++p)

#define FOREACH_NP(np) \
int zStart = max(k-1,0), zEnd = min(k+1,gwmm);\
int yStart = max(j-1,0), yEnd = min(j+1,ghmm);\
int xStart = max(i-1,0), xEnd = min(i+1,gdmm);\
int np, npEnd, gridZ, gridYZ, gridXYZ, gridZEnd, gridYZEnd, gridXYZEnd;\
int gridYStart = gw * yStart;\
int gridYEnd = gw * yEnd;\
for (gridZ = gwh * zStart, gridZEnd = gwh * zEnd; gridZ <= gridZEnd; gridZ += gwh)\
for (gridYZ = gridZ + gridYStart, gridYZEnd = gridZ + gridYEnd; gridYZ <= gridYZEnd; gridYZ += gw)\
for (gridXYZ = gridYZ + xStart, gridXYZEnd = gridYZ + xEnd; gridXYZ <= gridXYZEnd; ++gridXYZ)\
for (np = gridStarts[gridXYZ], npEnd = gridEnds[gridXYZ]; np < npEnd; ++np) \
if (p <= np)

void FluidSystem::updateDensities()
{
	FOREACH_P(p) {
		float pDensity = 0.f;
		float4 particlePosition = particlePositions[p];
		FOREACH_NP(np) {
			float4 d = particlePosition - particlePositions[np];
			float rSq = d.absSquared();
			if (rSq <= HH) {
				float kernelDiff = HH - rSq;
				float common = KERNEL(kernelDiff) * MASS;
				pDensity += common;
				particleDensities[np] += common;
			}
		}
		particleDensities[p] += pDensity;
	}
}



void FluidSystem::updateForces()
{
	float4 gravCommon = MASS * GRAV_CONST * gravDirection[0];
	FOREACH_P(p) {
		
		float particleDensity = particleDensities[p];
		float4 pForce = particleDensity * gravCommon;
		float4 particlePosition = particlePositions[p];
		float particleDensityRecip = MASS / particleDensity;
		
		FOREACH_NP(np) {
			float4 d = particlePosition - particlePositions[np];
			float rSq = d.absSquared();
			if (rSq <= HH) {
				float r = sqrtf(rSq);
				float neighborDensity = particleDensities[np];
				float neighborDensityRecip = MASS / neighborDensity;
				float kernelDiff = HH - rSq;
				float kernelDiff2 = H - r;
				float4 npForce;
				
				// Compute the pressure force.
				float4 common = HALF_K
					* (particleDensity + neighborDensity - REST_DENSITY_TIMES_2)
					* PRESSURE_GRADIENT(d, r, kernelDiff2);
				pForce  -= neighborDensityRecip * common;
				npForce += particleDensityRecip * common;
				
				// Compute the viscosity force.
				common = MU * (particleVelocities[np] - particleVelocities[p])
					* VISCOSITY_LAPLACIAN(kernelDiff2);
				pForce  += neighborDensityRecip * common;
				npForce -= particleDensityRecip * common;
				
				// Compute the gradient of the color field.
				common = GRADIENT(d, kernelDiff);
				float4 particleColorGradient  = neighborDensityRecip * common;
				float4 neighbourColorGradient = particleDensityRecip * common;
				
				// Compute the laplacian of the color field.
				float value = SIGMA * LAPLACIAN(kernelDiff);
				
				float particleColorGradLen = particleColorGradient.abs();
				if (particleColorGradLen > 0.001f) {
					pForce -= neighborDensityRecip * value
					* (particleColorGradient / particleColorGradLen);
				}
				
				float neighborColorGradLen = neighbourColorGradient.abs();
				if (neighborColorGradLen > 0.001f) {
					npForce -= particleDensityRecip * value * (neighbourColorGradient / neighborColorGradLen);
				}
				particleForces[np] += npForce;
				
			}
		}
		particleForces[p] += pForce;
	}
	
}

void FluidSystem::updateParticles()
{
	FOREACH_PARTICLE(p) {
		float4 vel = particleVelocities[p];
		float4 pos = particlePositions[p];
		vel += TIME_STEP * (particleForces[p] / particleDensities[p] - (POINT_DAMPING * vel ) / MASS); // step vel
		pos += TIME_STEP * vel; // step pos
		
		// boundaries
		if (pos[0] < 0.f) {
			pos[0] = 0.f;
			vel[0] *= -COLLISION_ELASTICITY;
		} else if (pos[0] > WIDTH) {
			pos[0] = WIDTH;
			vel[0] *= -COLLISION_ELASTICITY;
		}
		
		if (pos[1] < 0.f) {
			pos[1] = 0.f;
			vel[1] *= -COLLISION_ELASTICITY;
		} else if (pos[1] > HEIGHT) {
			pos[1] = HEIGHT;
			vel[1] *= -COLLISION_ELASTICITY;
		}
		
		if (pos[2] < 0.f) {
			pos[2] = 0.f;
			vel[2] *= -COLLISION_ELASTICITY;
		} else if (pos[2] > DEPTH) {
			pos[2] = DEPTH;
			vel[2] *= -COLLISION_ELASTICITY;
		}
		
		particlePositions[p] = pos; // write back pos
		particleVelocities[p] = vel; // write back vel
		
	}
}

void FluidSystem::draw()
{
#ifdef OPENCL_SPH
	// Calculate
	CL_ENQUEUE_WRITE_BUFFER(queue, cl_gravDirection, sizeof(float4), gravDirection);
	
	CL_ENQUEUE_ND_RANGE_KERNEL_AUTO_LOCAL(queue, updateDensityKernel, 1, PARTICLE_COUNT);
	CL_ENQUEUE_ND_RANGE_KERNEL_AUTO_LOCAL(queue, updateForceKernel, 1, PARTICLE_COUNT);
	CL_ENQUEUE_ND_RANGE_KERNEL_AUTO_LOCAL(queue, updateParticleKernel, 1, PARTICLE_COUNT);
#else
	updateDensities();
	updateForces();
	updateParticles();
#endif
	
	// Draw
	if (!sphereInited) {
		sphereId = glGenLists(1);
		glNewList(sphereId, GL_COMPILE);
		glutSolidSphere(0.3 * SCALE, 5, 5);
		glEndList();
		sphereInited = true;
	}
	if (drawParticlesMode) {
		FOREACH_PARTICLE(p) drawParticle(p);
		if (++drawParticlesTimer > DRAW_PARTICLES_TIME_LIMIT) drawParticlesMode = false;
	} else {
		verts.clear();
		norms.clear();
		faces.clear();
		
		generateIsoSurface(verts, norms, faces);
		drawFaces(verts, norms, faces);
		
	}
#ifdef OPENCL_SPH
	CL_ENQUEUE_READ_BUFFER(queue, cl_particlePositions, sizeof(float4) * PARTICLE_COUNT, particlePositions);
	CL_ENQUEUE_READ_BUFFER(queue, cl_particleVelocities, sizeof(float4) * PARTICLE_COUNT, particleVelocities);
	CL_ENQUEUE_READ_BUFFER(queue, cl_particleDensities, sizeof(float) * PARTICLE_COUNT, particleDensities);
#endif
	updateGrid();
	
#ifdef OPENCL_SPH
	CL_ENQUEUE_WRITE_BUFFER(queue, cl_gridStarts, sizeof(int) * GRID_SIZE, gridStarts);
	CL_ENQUEUE_WRITE_BUFFER(queue, cl_gridEnds, sizeof(int) * GRID_SIZE, gridEnds);
	CL_ENQUEUE_WRITE_BUFFER(queue, cl_particlePositions, sizeof(float4) * PARTICLE_COUNT, particlePositions);
	CL_ENQUEUE_WRITE_BUFFER(queue, cl_particleVelocities, sizeof(float4) * PARTICLE_COUNT, particleVelocities);
#else
	FOREACH_PARTICLE(p) { particleDensities[p] = 0.f; particleForces[p] = 0.f; }
	
#endif
}


#ifdef OPENCL_SPH
#define GS prevGridStarts
#define PP gridMakeParticlePositions
#else
#define GS gridStarts
#define PP particlePositions
#endif

#define ADD_DENSITY_AT_POINT_BY_PARTICLE(var, point, p) \
{ float t = max(0.f, HH - (PP[p] - point).absSquared());\
var += particleDensities[p++] * CUBE(t); }

#define ADD_DENSITY_2 \
{ gridXYZ = gridYZ + xStart; \
p0 = GS[gridXYZ]; p1 = GS[++gridXYZ]; p2 = GS[++gridXYZ]; p3 = GS[++gridXYZ]; \
while (p0 < p1) ADD_DENSITY_AT_POINT_BY_PARTICLE(density, point, p0);\
while (p1 < p2) ADD_DENSITY_AT_POINT_BY_PARTICLE(density, point, p1);\
if (gridXYZ <= gridYZ + xEnd) while (p2 < p3) ADD_DENSITY_AT_POINT_BY_PARTICLE(density, point, p2); }
	
#define ADD_DENSITY \
{ gridYZ = gridZ + gridYStart; gridYZEnd = gridZ + gridYEnd; \
ADD_DENSITY_2; gridYZ += gw; \
if (density <= mct && gridYZ <= gridYZEnd) ADD_DENSITY_2; gridYZ += gw; \
if (density <= mct && gridYZ <= gridYZEnd) ADD_DENSITY_2; }

#define SET_DENSITY(var, p) \
{ float4 point = p;\
float4 pointIndex = (point / H).floor();\
float density = 0.f;\
float mct = MC_THRESHOLD;\
\
static float4 one = float4(1.f);\
float4 xyzStart = max(pointIndex - one, float4());\
float4 xyzEnd = min(pointIndex + one, float4(GRID_WIDTH, GRID_HEIGHT, GRID_DEPTH) - one);\
\
int xStart = xyzStart[0]; int xEnd = xyzEnd[0];\
int yStart = xyzStart[1]; int yEnd = xyzEnd[1];\
int zStart = xyzStart[2]; int zEnd = xyzEnd[2];\
\
int p0, p1, p2, p3, gridZ, gridYZ, gridYZEnd, gridXYZ, gridZEnd;\
int gridYStart = GRID_WIDTH * yStart;\
int gridYEnd = GRID_WIDTH * yEnd;\
int gw = GRID_WIDTH;\
int gwh = GRID_WIDTH_HEIGHT;\
\
gridZ = gwh * zStart; gridZEnd = gwh * zEnd;\
ADD_DENSITY; gridZ += gwh;\
if (density <= mct && gridZ <= gridZEnd) ADD_DENSITY; gridZ += gwh;\
if (density <= mct && gridZ <= gridZEnd) ADD_DENSITY;\
\
var = density; }

#define DO_8x(c) c;c;c;c;c;c;c;c;

inline void FluidSystem::generateIsoSurface(float4v &verts, float4v &norms, mini_vec<Face> &faces)
{
	static const int m = MC_MARGIN;
	
	static int stepsX = int(WIDTH / MC_RES) ;
	static int stepsY = int(HEIGHT / MC_RES) ;
	static int stepsZ = int(DEPTH / MC_RES) ;
	static int stepsYPP = stepsY + m + m + stepsX;
	static int stepsZPP = stepsZ + m + m;
	static int stepsYZPP = stepsYPP * stepsZPP;
	
	static float4 cornerOffsetVecs[8];
	static const int cornerOffsets[8][3] =  {
		{0,0,1}, {1,0,1}, {1,0,0}, {0,0,0},
		{0,1,1}, {1,1,1}, {1,1,0}, {0,1,0} };
	
	int i=0;
	DO_8x(
		cornerOffsetVecs[i] = float4(cornerOffsets[i][0],
									 cornerOffsets[i][1],
									 cornerOffsets[i][2]);
		++i;
	)
	float *isoCache = new float[2 * stepsYZPP + 10];
	CubeEdges *cubeEdges = new CubeEdges[2 * stepsYZPP + 10];
	
	CubeEdges emptyCube;
	
	int currentCubeEdgesLayer = 0;
	int previousCubeEdgesLayer = 1;
	
	bool firstRun = true;
	int nextIsoCacheLayer = 0;
	
	static int negM = -m;
	static int xEnd = stepsX + m;
	static int yEnd = stepsY + m;
	static int zEnd = stepsZ + m;
	
	for (int x=negM; x<xEnd; ++x) {
		int isoCacheLayers[2] = { (nextIsoCacheLayer + 1) & 1, nextIsoCacheLayer };
		nextIsoCacheLayer = isoCacheLayers[0];
		
		for (int r=0; r<1+firstRun; ++r) {
			int cacheBase = isoCacheLayers[1-r] * stepsYZPP;
			int xpp = x+1-r;
			for (int y=negM; y<=yEnd; ++y) {
				int cacheBaseY = cacheBase + (y+m) * stepsZPP;
				for (int z=negM; z<=zEnd; ++z) {
					SET_DENSITY(isoCache[cacheBaseY + z + m], float4(xpp, y, z) * MC_RES);
				}
			}
		}
		
		int cl = currentCubeEdgesLayer * stepsYZPP;
		int pl = previousCubeEdgesLayer * stepsYZPP;
		
		for (int y=negM; y<yEnd; ++y) {
			int ypp = y+m;
			for (int z=negM; z<zEnd; ++z) {
				Cube cell;
				int i=0, cacheBase, cacheBaseY;
				float4 xyz(x,y,z);
				int zpp = z+m;
				DO_8x(
					cell.p[i] = SCALE * ((xyz + cornerOffsetVecs[i]) * MC_RES + drawingOffset);
					cacheBase = isoCacheLayers[cornerOffsets[i][0]] * stepsYZPP;
					cacheBaseY = cacheBase + (cornerOffsets[i][1] + ypp) * stepsZPP;
					cell.val[i] = isoCache[cacheBaseY + zpp + cornerOffsets[i][2]];
					++i;
				)
				int yl = ypp*stepsYPP;
				CubeEdges &c0 = z>negM ? cubeEdges[cl + yl + zpp - 1] : emptyCube;
				CubeEdges &c1 = y>negM ? cubeEdges[cl + (ypp - 1) * stepsYPP + zpp] : emptyCube;
				CubeEdges &c2 = x>negM ? cubeEdges[pl + yl + zpp] : emptyCube;
				cubeEdges[cl + ypp * stepsYPP + zpp] =
					polygonize(cell, MC_POLY_THRESHOLD, c0, c1, c2, verts, faces);
				
			}
		}
		
		firstRun = false;
		swap(currentCubeEdgesLayer, previousCubeEdgesLayer);
	}
	delete [] isoCache;
	delete [] cubeEdges;

	calculateNorms(verts, norms, faces);
	
}

void FluidSystem::setDrawParticleMode(bool mode)
{
	drawParticlesMode = mode;
	if (drawParticlesMode) drawParticlesTimer = 0;
}

bool FluidSystem::getDrawParticleMode()
{
	return drawParticlesMode;
}

