#ifndef SHADER_HELPERS_H
#define SHADER_HELPERS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _WIN32
#include "GL/freeglut.h"
#elif __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "bitmap_image.hpp"

inline char *textFileRead(const char *fn)
{
	FILE *fp;
	char *content = NULL;
	size_t count=0;
	
	if (fn != NULL) {
		fp = fopen(fn,"rt");
		
		if (fp != NULL) {
			fseek(fp, 0, SEEK_END);
			count = ftell(fp);
			rewind(fp);
			
			if (count > 0) {
				content = (char *)malloc(sizeof(char) * (count+1));
				count = fread(content,sizeof(char),count,fp);
				content[count] = '\0';
			}
			fclose(fp);
		}
	}
	return content;
}


GLuint makeGLProgram(const char* vertFilename, const char* fragFilename)
{
	char *vs = NULL, *fs = NULL;
	
	GLuint glProgram;
	GLuint vertShader = glCreateShader(GL_VERTEX_SHADER);
	GLuint fragShader = glCreateShader(GL_FRAGMENT_SHADER);
	
	vs = textFileRead(vertFilename);
	fs = textFileRead(fragFilename);
	
	const char * ff = fs;
	const char * vv = vs;
	
	glShaderSource(vertShader, 1, &vv, NULL);
	glShaderSource(fragShader, 1, &ff, NULL);
	
	free(vs);free(fs);
	
	glCompileShader(vertShader);
	glCompileShader(fragShader);
	
	int status;
	glGetShaderiv(vertShader, GL_COMPILE_STATUS, &status);
	if (status != GL_TRUE) {
		GLint maxLength = 0;
		glGetShaderiv(vertShader, GL_INFO_LOG_LENGTH, &maxLength);
		GLchar errorLog[maxLength];
		glGetShaderInfoLog(vertShader, maxLength, &maxLength, errorLog);
		printf("Vertex shader ""%s"" failed to compile. \n%s\n", vertFilename, errorLog);
	}
	glGetShaderiv(fragShader, GL_COMPILE_STATUS, &status);
	if (status != GL_TRUE) {
		GLint maxLength = 0;
		glGetShaderiv(fragShader, GL_INFO_LOG_LENGTH, &maxLength);
		GLchar errorLog[maxLength];
		glGetShaderInfoLog(fragShader, maxLength, &maxLength, errorLog);
		printf("Fragment shader ""%s"" failed to compile. \n%s\n", fragFilename, &errorLog[0]);
	}
	
	glProgram = glCreateProgram();
	glAttachShader(glProgram,fragShader);
	glAttachShader(glProgram,vertShader);
	
	glLinkProgram(glProgram);
	glGetProgramiv(glProgram, GL_LINK_STATUS, &status);
	if (status != GL_TRUE) {
		printf("Failed to link ""%s"" and ""%s"".\n", vertFilename, fragFilename);
	}
	
	glDeleteShader(vertShader);
	glDeleteShader(fragShader);
	
	return glProgram;
}

GLuint makeCubemap(const char *rightFilename,
				   const char *leftFilename,
				   const char *topFilename,
				   const char *bottomFilename,
				   const char *backFilename,
				   const char *frontFilename,
				   bool bgr=false)
{
	GLuint textureID;
	glGenTextures(1, &textureID);
	
	bitmap_image face;
	glBindTexture(GL_TEXTURE_CUBE_MAP, textureID);
	
#define LOAD_FACE(filename, index) \
face = bitmap_image(filename); \
glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_X + index, 0, GL_RGB,\
face.width(), face.height(), 0, bgr?GL_BGR:GL_RGB, GL_UNSIGNED_BYTE, face.data());
	LOAD_FACE(rightFilename, 0);
	LOAD_FACE(leftFilename, 1);
	LOAD_FACE(topFilename, 2);
	LOAD_FACE(bottomFilename, 3);
	LOAD_FACE(backFilename, 4);
	LOAD_FACE(frontFilename, 5);
#undef LOAD_FACE
	
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
	glBindTexture(GL_TEXTURE_CUBE_MAP, 0);
	
	return textureID;
}

GLuint makeSkyboxDisplayList(GLfloat size)
{
	GLfloat skyboxVertices[] = {
		-1.0f, 1.0f, -1.0f,
		-1.0f, -1.0f, -1.0f,
		1.0f, -1.0f, -1.0f,
		1.0f, -1.0f, -1.0f,
		1.0f, 1.0f, -1.0f,
		-1.0f, 1.0f, -1.0f,
		
		-1.0f, -1.0f, 1.0f,
		-1.0f, -1.0f, -1.0f,
		-1.0f, 1.0f, -1.0f,
		-1.0f, 1.0f, -1.0f,
		-1.0f, 1.0f, 1.0f,
		-1.0f, -1.0f, 1.0f,
		
		1.0f, -1.0f, -1.0f,
		1.0f, -1.0f, 1.0f,
		1.0f, 1.0f, 1.0f,
		1.0f, 1.0f, 1.0f,
		1.0f, 1.0f, -1.0f,
		1.0f, -1.0f, -1.0f,
		
		-1.0f, -1.0f, 1.0f,
		-1.0f, 1.0f, 1.0f,
		1.0f, 1.0f, 1.0f,
		1.0f, 1.0f, 1.0f,
		1.0f, -1.0f, 1.0f,
		-1.0f, -1.0f, 1.0f,
		
		-1.0f, 1.0f, -1.0f,
		1.0f, 1.0f, -1.0f,
		1.0f, 1.0f, 1.0f,
		1.0f, 1.0f, 1.0f,
		-1.0f, 1.0f, 1.0f,
		-1.0f, 1.0f, -1.0f,
		
		-1.0f, -1.0f, -1.0f,
		-1.0f, -1.0f, 1.0f,
		1.0f, -1.0f, -1.0f,
		1.0f, -1.0f, -1.0f,
		-1.0f, -1.0f, 1.0f,
		1.0f, -1.0f, 1.0f
	};
	GLuint skyboxDispList = glGenLists(1);
	glNewList(skyboxDispList, GL_COMPILE);

	glBegin(GL_TRIANGLES);
	for (int i=0; i<36; ++i) {
		glVertex3f(skyboxVertices[i*3] * size, skyboxVertices[i*3+1] * size, skyboxVertices[i*3+2] * size);
	}
	glEnd();
	
	glEndList();
	
	return skyboxDispList;
}

#define GL_UNIFORM_3FV(program, string, value) \
glUniform3fv(glGetUniformLocation(program, string), 1, value);

#define GL_UNIFORM_1I(program, string, value) \
glUniform1i(glGetUniformLocation(program, string), value);

#define GL_UNIFORM_MAT_4FV(program, string, value) \
glUniformMatrix4fv(glGetUniformLocation(program, string), 1, GL_FALSE, value);

#define GL_UNIFORM_MAT_3FV(program, string, value) \
glUniformMatrix3fv(glGetUniformLocation(program, string), 1, GL_FALSE, value);


#endif
