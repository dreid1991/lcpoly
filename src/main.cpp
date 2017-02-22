#include "time.h"
#include <iostream>
#include <math.h>
#include <cstdlib>
#include <stdio.h>
#include "Sim.hpp"

int main() {
	Sim *mySim = new Sim();
	mySim->MCSim();
	
    free(mySim);
	puts("SIMULATION COMPLETE");
	
	
	return 0;
}
