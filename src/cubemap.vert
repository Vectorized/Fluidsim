varying vec3 refracted, reflected;
varying float schlick_R;
uniform vec3 cameraCenter;

// precalculated for water with refr index of 1.33
#define RATIO 0.7500001875
#define R_0 0.0204081282799

uniform mat3 reverseMapMat;
uniform mat4 proj;
uniform mat3 normalMat;
uniform mat4 modelView;
uniform mat4 modelViewProj;

void main()
{
	vec3 normal, position, I;
	float mc, mc2;
	
	normal = normalize(normalMat * gl_Normal);
	position = vec3(modelView * gl_Vertex);

	I = normalize(position - cameraCenter);

	mc = 1.0 - abs(dot(I, normal));
	mc2 = mc * mc;
	schlick_R = clamp(R_0 + (1.0 - R_0) * mc2 * mc2 * mc, 0.0, 1.0);
	
	refracted = reverseMapMat * normalize(refract(I, normal, RATIO));
	reflected = reverseMapMat * normalize(reflect(I, normal));
	
	gl_Position = modelViewProj * gl_ModelViewMatrix * gl_Vertex;
}