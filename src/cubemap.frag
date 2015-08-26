varying vec3 refracted, reflected;
varying float schlick_R;

uniform samplerCube skybox;

void main()
{
	gl_FragColor = (1.0 - schlick_R) * textureCube(skybox, refracted) +
					schlick_R * textureCube(skybox, reflected);
}