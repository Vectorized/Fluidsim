#ifndef FPS_H
#define FPS_H

#include <time.h>

struct FPS
{
	int windowWidth, windowHeight;
	int frameCount;
	int textX, textY;
	float fps;
	int currentTime, previousTime;
	void *fontStyle;
	
	FPS(int windowWidth, int windowHeight, int textX, int textY, void *fontStyle=GLUT_BITMAP_HELVETICA_12) {
		this->windowWidth = windowWidth;
		this->windowHeight = windowHeight;
		this->textX = textX;
		this->textY = textY;
		frameCount = currentTime = previousTime = 0;
		fps = 0.f;
		this->fontStyle = fontStyle;
	}
	
	void updateWindowSize(int windowWidth, int windowHeight) {
		this->windowWidth = windowWidth;
		this->windowHeight = windowHeight;
	}

	
	void drawBitmapText(char *string, float x, float y, float z) {
		char *c;
		glRasterPos3f(x, y, z);
		for (c=string; *c != '\0'; c++)
			glutBitmapCharacter(fontStyle, *c);
	}
	
	void draw() {
		frameCount++;
		currentTime = glutGet(GLUT_ELAPSED_TIME);
		int timeInterval = currentTime - previousTime;
		
		if(timeInterval > 1000) {
			fps = frameCount / (timeInterval / 1000.0f);
			previousTime = currentTime;
			frameCount = 0;
		}
		
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glOrtho( 0, windowWidth, 0, windowHeight, -1, 1 );
		
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		
		glDisable( GL_LIGHTING );
		glColor3ub(255,255,255);
		char fpsStr[1024];
		snprintf(fpsStr, sizeof(fpsStr), "Frame rate: %.2f", fps);
		drawBitmapText(fpsStr, textX, textY, 0);
		glEnable( GL_LIGHTING );
		
	}
};


#endif
