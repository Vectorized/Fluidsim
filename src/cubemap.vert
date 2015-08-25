varying vec3 normal, position;

uniform mat4 proj;
uniform mat3 normalMat;
uniform mat4 modelView;
uniform mat4 modelViewProj;

void main()
{
	normal = normalize(normalMat * gl_Normal);
	position = vec3(modelView * gl_Vertex);
	gl_Position = modelViewProj * gl_ModelViewMatrix * gl_Vertex;
}