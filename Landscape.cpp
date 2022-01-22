#include "Landscape.h"

Landscape::Landscape() {
	pres = 0;
	suitable = 0;
	local_K = 0.0; 
	local_K0 = NULL;  
	eps = 0.0;
	//theta = NULL; 
}

//------------------------------------------------------------------
Landscape::~Landscape() {
	if (local_K0 != NULL) { delete local_K0; local_K0 = NULL; }
	//delete theta; theta = NULL; 
}
