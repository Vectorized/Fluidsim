#ifndef GnV_FluidSim_exit_pass_h
#define GnV_FluidSim_exit_pass_h

#include <string>

#define KIOSK_MODE

class ExitPass
{
private:
	bool active;
	int currentIndex;
	string pass;
	
public:
	ExitPass(string pass)
	{
		active = false;
		currentIndex = 0;
		this->pass = pass;
	}
	
	void activate()
	{
#ifdef KIOSK_MODE
		active = 1;
#else
		exit(0);
#endif
	}
	
	bool isActive()
	{
		return active;
	}
	
	void enterChar(unsigned char c)
	{
		if (c == 27) {
			activate();
			return;
		}
		
		if (active && c == pass.c_str()[currentIndex]) {
			if (++currentIndex == pass.size()) exit(0);
		} else {
			active = 0;
			currentIndex = 0;
		}
	}
	
};


#endif
