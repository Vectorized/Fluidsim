varying vec3 texCoords;

uniform mat4 viewProj;

void main()
{
	texCoords = vec3(gl_Vertex);
	gl_Position = viewProj * gl_Vertex;
}