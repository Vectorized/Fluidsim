#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <vector>

#ifdef _WIN32
#include "glew.h"
#include "wglew.h"
#include "GL/freeglut.h"
#elif __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "vecmath.h"
#include "camera.h"
#include "fluidSystem.h"
#include "shader_helpers.h"
#include "fps.h"
#include "showcase.h"
#include "exit_pass.h"

using namespace std;

// Globals here.
const bool showText = true;
const int windowHeight = 600;
const int windowWidth = 600;

GLuint skyboxProgram,
skyboxTexture,
skyboxDisplayListLarger,
skyboxDisplayListSmaller,
cubemapProgram;


FluidSystem f(18, 18, 18);

// This is the camera
Camera camera;

FPS fps(windowWidth, windowHeight, 10, 10);

Instructions instructions = Instructions(
	"Instructions:\n"
	"Use the mouse to drag and rotate the water's box.\n"
	"Use the arrow keys to rotate the view.\n"
	"Press 'm' to toggle render mode.",
	windowWidth, windowHeight, 10, windowHeight - 20, 18);

AutoRotator autoRotator = AutoRotator(&camera);

ExitPass exitPass("pass");

//-------------------------------------------------------------------
// These are state variables for the UI
bool g_mousePressed = false;

// Declarations of functions whose implementations occur later.
void keyboardFunc( unsigned char key, int x, int y);
void specialFunc( int key, int x, int y );
void mouseFunc(int button, int state, int x, int y);
void motionFunc(int x, int y);
void reshapeFunc(int w, int h);
void drawScene(void);
void initRendering();

// This function is called whenever a "Normal" key press is
// received.
void keyboardFunc(unsigned char key, int x, int y)
{
	if (exitPass.isActive()) {
		exitPass.enterChar(key);
		return;
	}
	
	if (key == 27) {
		exitPass.activate();
	} else if (key == ' ') {
		Matrix4f eye = Matrix4f::identity();
		camera.SetRotation( eye );
		camera.SetCenter( Vector3f::ZERO );
	} else if (key == 'm') {
		if (f.getDrawParticleMode()) {
			f.setDrawParticleMode(false);
			cout << "Particle rendering off... using Marching Cubes to draw liquid isosurface." << endl;
		} else {
			f.setDrawParticleMode(true);
			cout << "Particle rendering on... drawing individual liquid SPH particles." << endl;
		}
	} else if (key == 'f') {
		glutFullScreen();
	} else {
		cout << "Unhandled key press " << key << "." << endl;
	}
	glutPostRedisplay();
}


enum { KEY_LEFT = 100, KEY_UP = 101, KEY_RIGHT = 102, KEY_DOWN = 103 };

// This function is called whenever a "Special" key press is
// received.  Right now, it's handling the arrow keys.
void specialFunc(int key, int x, int y)
{
	float step = 0.025f;
	if (key == KEY_LEFT) {
		camera.ViewArcBallRotateAboutYDown(step);
	} else if (key == KEY_UP) {
		camera.ViewArcBallRotateAboutXDown(step);
	} else if (key == KEY_RIGHT) {
		camera.ViewArcBallRotateAboutYDown(-step);
	} else if (key == KEY_DOWN) {
		camera.ViewArcBallRotateAboutXDown(-step);
	}
}

void specialUpFunc(int key, int x, int y)
{
	if (key == KEY_LEFT || key == KEY_RIGHT) camera.ViewArcBallRotateAboutYRelease();
	else if (key == KEY_UP || key == KEY_DOWN) camera.ViewArcBallRotateAboutXRelease();
}

//  Called when mouse button is pressed.
void mouseFunc(int button, int state, int x, int y)
{
	if (state == GLUT_DOWN) {
		g_mousePressed = true;
		autoRotator.resetIdle();
		autoRotator.idleCounterOn = false;
#define LEFT_CLICK_ONLY
#ifdef LEFT_CLICK_ONLY
		camera.MouseClick(Camera::LEFT, x, y);
#else
		switch (button) {
			case GLUT_LEFT_BUTTON:
				camera.MouseClick(Camera::LEFT, x, y);
				break;
			case GLUT_MIDDLE_BUTTON:
				camera.MouseClick(Camera::MIDDLE, x, y);
				break;
			case GLUT_RIGHT_BUTTON:
				camera.MouseClick(Camera::RIGHT, x,y);
			default:
				break;
		}
#endif
	} else {
		autoRotator.idleCounterOn = true;
		camera.MouseRelease(x,y);
		g_mousePressed = false;
	}
	glutPostRedisplay();
}

// Called when mouse is moved while button pressed.
void motionFunc(int x, int y)
{
	camera.MouseDrag(x,y);
	glutPostRedisplay();
}

// Called when the window is resized
// w, h - width and height of the window in pixels.
void reshapeFunc(int w, int h)
{
	autoRotator.resetIdle();
	
	camera.SetDimensions(w,h);
	camera.SetViewport(0,0,w,h);
	camera.ApplyViewport();
	camera.SetPerspective(50);
}

// Initialize OpenGL's rendering modes
void initRendering()
{
	glEnable(GL_DEPTH_TEST);

	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	
	glClearColor(0,0,0,1);
	
	skyboxProgram = makeGLProgram("skybox.vert", "skybox.frag");
	skyboxTexture = makeCubemap("right.bmp",
								"left.bmp",
								"top.bmp",
								"bottom.bmp",
								"back.bmp",
								"front.bmp",
								true);
	
	skyboxDisplayListLarger  = makeSkyboxDisplayList(100);
	skyboxDisplayListSmaller = makeSkyboxDisplayList(50);
	
	cubemapProgram = makeGLProgram("cubemap.vert", "cubemap.frag");
	
	glActiveTexture(GL_TEXTURE0);
	
	GL_UNIFORM_1I(cubemapProgram, "skybox", 0);
	GL_UNIFORM_1I(skyboxProgram, "skybox", 0);

	glBindTexture(GL_TEXTURE_CUBE_MAP, skyboxTexture);
}


// This function is responsible for displaying the object.
void drawScene(void)
{
	// Clear the rendering window
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
	glMatrixMode( GL_MODELVIEW );
	glLoadIdentity();
	
	// Light color (RGBA)
	GLfloat Lt0diff[] = {1.0,1.0,1.0,1.0};
	GLfloat Lt0pos[] = {3.0,3.0,5.0,1.0};
	glLightfv(GL_LIGHT0, GL_DIFFUSE, Lt0diff);
	glLightfv(GL_LIGHT0, GL_POSITION, Lt0pos);
	
	M4 view  = camera.viewMatrix();
	M4 model = camera.modelMatrix();
	M4 proj  = camera.projectionMatrix();
	M4 modelView = view * model;
	M4 viewProj = proj * view;
	M4 modelViewProj = proj * view * model;
	M3 normalMat = modelView.inverse().transposed().getSubmatrix3x3(0, 0);
	M3 reverseMapMat = view.inverse().getSubmatrix3x3(0, 0);
	
	// DRAW MODEL START ======================
	
	glUseProgram(cubemapProgram);
	
	V3 cameraCenter = camera.GetCenter();

	GL_UNIFORM_MAT_4FV(cubemapProgram, "proj", proj);
	GL_UNIFORM_MAT_3FV(cubemapProgram, "normalMat", normalMat);
	GL_UNIFORM_MAT_4FV(cubemapProgram, "modelView", modelView);
	GL_UNIFORM_MAT_4FV(cubemapProgram, "modelViewProj", modelViewProj);
	GL_UNIFORM_MAT_3FV(cubemapProgram, "reverseMapMat", reverseMapMat);
	GL_UNIFORM_3FV(cubemapProgram, "cameraCenter", cameraCenter);

	glCallList(skyboxDisplayListLarger);
	
	f.setGravDirection(float4(-model[1], -model[5], -model[9]).normalized());
	f.draw();
	
	// DRAW MODEL END ========================
	
	// DRAW SKYBOX START =====================

	glDepthMask(GL_FALSE);
	
	glUseProgram(skyboxProgram);

	view(0,3) = view(1,3) = view(2,3) = 0.f;
	GL_UNIFORM_MAT_4FV(skyboxProgram, "viewProj", viewProj);
	glCallList(skyboxDisplayListSmaller);
	glDepthMask(GL_TRUE);

	glUseProgram(0);
 
	// DRAW SKYBOX END =======================
	
	autoRotator.stepAutoRotator();
	camera.ViewArcBallRotateStep();
	if (showText) {
		fps.draw();
		instructions.draw();
	}
	
	// Dump the image to the screen.
	glutSwapBuffers();
}


void timerFunc(int t)
{
	glutPostRedisplay();
	glutTimerFunc(t, &timerFunc, t);
}

// Main routine.
// Set up OpenGL, define the callbacks and start the main loop
int main( int argc, char* argv[] )
{

	glutInit( &argc, argv );
	
	// We're going to animate it, so double buffer
	glutInitDisplayMode( GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH );
	
	// Initial parameters for window position and size
	glutInitWindowPosition( 60, 60 );
	glutInitWindowSize( windowWidth, windowHeight );
	
	camera.SetDimensions( windowWidth, windowHeight );
	
	camera.SetDistance( 10 );
	camera.SetCenter( Vector3f::ZERO );
	
	glutCreateWindow("50.017 Fluid Sim");
	
	// Initialize OpenGL parameters.
#ifdef _WIN32
	glewInit();
#endif
	initRendering();
	
	// Set up callback functions for key presses
	glutKeyboardFunc(keyboardFunc);   // Handles "normal" ascii symbols
	glutSpecialFunc(specialFunc);     // Handles "special" keyboard keys on pressed down
	glutSpecialUpFunc(specialUpFunc); // Handles "special" keyboard keys on release
	glutMouseFunc(mouseFunc);
	glutMotionFunc(motionFunc);
	glutReshapeFunc( reshapeFunc );
	glutDisplayFunc( drawScene );
	glutTimerFunc(15, timerFunc, 15);
	
	// Start the main loop.  glutMainLoop never returns.
	glutMainLoop();
		
	return 0;	// This line is never reached.
}
