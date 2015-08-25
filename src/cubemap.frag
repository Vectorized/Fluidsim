varying vec3 normal, position;

uniform vec3 cameraCenter;
uniform samplerCube skybox;

// precalculated for water with refr index of 1.33
#define RATIO 0.7500001875
#define R_0 0.0204081282799

uniform mat3 reverseMapMat;


void main()
{
	float mc;
	float schlick_R;
	vec3 I;

	I = normalize(position - cameraCenter);
	
	mc = 1.0 - abs(dot(I, normal));
	schlick_R = clamp(R_0 + (1.0 - R_0) * mc * mc * mc * mc * mc, 0.0, 1.0);
	
	gl_FragColor = (1.0 - schlick_R) * textureCube(skybox, reverseMapMat * normalize(refract(I, normal, RATIO))) +
					schlick_R * textureCube(skybox, reverseMapMat * normalize(reflect(I, normal)));
}