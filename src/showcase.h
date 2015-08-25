#ifndef GnV_FluidSim_showcase_h
#define GnV_FluidSim_showcase_h

#include "camera.h"
#include "helpers.h"

struct AutoRotator
{
	int idleTime;
	Camera *c;
	M4 rotateMat;
	bool idleCounterOn;
	AutoRotator(Camera *c) {
		idleCounterOn = true;
		idleTime = 0;
		this->c = c;
		rotateMat = M4::rotateX(0.008)*M4::rotateZ(0.0038);
	}
	void resetIdle() {
		idleTime = 0;
	}
	void stepAutoRotator() {
		if (idleTime < 118) {
			idleTime += idleCounterOn;
		} else {
			c->SetRotation(rotateMat * c->GetRotation());
		}
	}
};

struct Instructions
{
	int windowWidth, windowHeight;
	int textX, textY, lineHeight;
	GLvoid *fontStyle;
	string text;
	
	Instructions(int windowWidth, int windowHeight, int textX, int textY, int lineHeight, GLvoid *fontStyle=GLUT_BITMAP_HELVETICA_12) {
		this->windowWidth = windowWidth;
		this->windowHeight = windowHeight;
		this->textX = textX;
		this->textY = textY;
		this->lineHeight = lineHeight;
		this->fontStyle = fontStyle;
		this->text = "Instructions:\n"
					 "Use the mouse to drag and rotate the water's box.\n"
					 "Use the arrow keys to rotate the view.\n"
					 "Press 'm' to toggle render mode.";
	}
	
	void updateWindowSize(int windowWidth, int windowHeight) {
		this->windowWidth = windowWidth;
		this->windowHeight = windowHeight;
	}

	
	void draw() {
		
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glOrtho( 0, windowWidth, 0, windowHeight, -1, 1 );
		
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		
		glDisable( GL_LIGHTING );
		
		glColor3ub(255,255,255);
		int i = 0;
		glRasterPos3f(textX, textY, 0);
		for (const char *c = text.c_str(); *c != '\0'; c++) {
			if (*c == '\n') {
				glRasterPos3f(textX, textY + --i * lineHeight, 0);
			} else {
				glutBitmapCharacter(fontStyle, *c);
			}
		}

		glEnable( GL_LIGHTING );
		
	}
};


#endif