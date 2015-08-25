varying vec3 texCoords;

uniform samplerCube skybox;

void main()
{
	gl_FragColor = textureCube(skybox, texCoords);
}

