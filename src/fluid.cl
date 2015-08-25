#define SQR(x)	((x)*(x))
#define CUBE(x)	((x)*(x)*(x))

#define GRID_INDEX(i,j,k) (GRID_WIDTH_HEIGHT * k + GRID_WIDTH * j + i)

#define KERNEL(kernelDiff) KERNEL_BASE * CUBE(kernelDiff)

#define GRADIENT(d, kernelDiff) (-GRADIENT_OR_LAPLACIAN_BASE * SQR(kernelDiff) * d)

#define LAPLACIAN(kernelDiff) (GRADIENT_OR_LAPLACIAN_BASE * kernelDiff * (7.0f * rSq - 3.0f * HH))

#define PRESSURE_GRADIENT(d, r, kernelDiff2) \
((r<0.001f) ? 0.f : (PRESSURE_GRADIDENT_BASE * SQR(kernelDiff2) * d / r))

#define VISCOSITY_LAPLACIAN(kernelDiff2) (VISCOSITY_LAPLACIAN_BASE * kernelDiff2)

__kernel void updateDensity(__global int *starts, __global int *ends, __global float4 *positions, __global float *densities)
{
	int p = get_global_id(0);
	
	float rSq, kernelDiff;
	float4 d;
	int gnp, npEnd;
	
	float particleDensity = 0.f;
	float4 particlePosition = positions[p];
	
	float4 particleIndex = particlePosition * H_RECIP;
	int i = particleIndex[0];
	int j = particleIndex[1];
	int k = particleIndex[2];
	
	int zStart = max(k-1, 0); int zEnd = min(k+1, GRID_WIDTH-1);
	int yStart = max(j-1, 0); int yEnd = min(j+1, GRID_HEIGHT-1);
	int xStart = max(i-1, 0); int xEnd = min(i+1, GRID_DEPTH-1);
	
	for (int z = zStart; z <= zEnd; ++z)
	for (int y = yStart; y <= yEnd; ++y)
	for (int x = xStart; x <= xEnd; ++x) {
		gnp = GRID_INDEX(x, y, z);
		npEnd = ends[gnp];
		for (int np=starts[gnp]; np<npEnd; ++np) {
			d = particlePosition - positions[np];
			rSq = dot(d,d);
			kernelDiff = HH - rSq;
			particleDensity += (rSq <= HH) ? KERNEL(kernelDiff) * MASS : 0.f;
		}
	}
	densities[p] = particleDensity;

}

__kernel void updateForce(__global int *starts, __global int *ends, __global float4 *positions, __global float4 *velocities, __global float4 *forces, __global float *densities, __constant float4 *gravDirection)
{
	int p = get_global_id(0);
	
	float4 particlePosition = positions[p];
	
	float4 particleIndex = particlePosition * H_RECIP;
	int i = particleIndex[0];
	int j = particleIndex[1];
	int k = particleIndex[2];
	
	int zStart = max(k-1, 0); int zEnd = min(k+1, GRID_WIDTH-1);
	int yStart = max(j-1, 0); int yEnd = min(j+1, GRID_HEIGHT-1);
	int xStart = max(i-1, 0); int xEnd = min(i+1, GRID_DEPTH-1);
	
	float4 d, particleColorGradient;
	float rSq, r, neighborDensity, massNeighborDensityRecip, particleColorGradLen, kernelDiff, kernelDiff2;
	int gnp, npStart, npEnd;
	
	float particleDensity = densities[p];
	float4 particleForce = particleDensity * MASS * GRAV_CONST * gravDirection[0];

	for (int z = zStart; z <= zEnd; ++z)
	for (int y = yStart; y <= yEnd; ++y)
	for (int x = xStart; x <= xEnd; ++x) {
		gnp = GRID_INDEX(x, y, z);
		npEnd = ends[gnp];
		for (int np=starts[gnp]; np<npEnd; ++np) {
			
			d = particlePosition - positions[np];
			rSq = dot(d,d);
			if (rSq <= HH) {
				
				r = sqrt(rSq);
				neighborDensity = densities[np];
				massNeighborDensityRecip = MASS / neighborDensity;
				
				kernelDiff = HH - rSq;
				kernelDiff2 = H - r;
				
				// Compute the pressure force.
				particleForce  += -massNeighborDensityRecip * HALF_K
				* (particleDensity + neighborDensity - REST_DENSITY_TIMES_2)
				* PRESSURE_GRADIENT(d, r, kernelDiff2);
				
				// Compute the viscosity force.
				particleForce  += MASS  * MU * (velocities[np] - velocities[p])
				* VISCOSITY_LAPLACIAN(kernelDiff2);
				
				// Compute the gradient of the color field.
				particleColorGradient  = massNeighborDensityRecip * GRADIENT(d, kernelDiff);
				
				particleColorGradLen = fast_length(particleColorGradient);
				
				particleForce += (particleColorGradLen > 0.001f) ?
				-SIGMA * massNeighborDensityRecip * LAPLACIAN(kernelDiff) * particleColorGradient / particleColorGradLen :
				0.f;
			}
			
		}
	}
	forces[p] = particleForce;

}


__kernel void updateParticle(__global float4 *positions, __global float4 *velocities, __global float4 *forces, __global float *densities)
{
	int p = get_global_id(0);
	
	float4 vel = velocities[p];
	float4 pos = positions[p];
	vel += TIME_STEP * (forces[p] / densities[p] - (POINT_DAMPING * vel ) * MASS_RECIP); // step vel
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
	
	positions[p] = pos; // write back pos
	velocities[p] = vel; // write back vel
}

