#include "defs.h"


class Disk {
	public:
        Vector3d r; //Contains x,y,z coordinates
        Vector3d rn; //Contains x,y,z coordinates without periodic boundary conditions for bond calculation purposes
        Vector3d u; //Vector tangent to the polymer 
        Vector3d f; //Vector normal to the disk (used for Gay-Berne calculations)
        Vector3d v; //Vector normal to u and f

        /* Using this to store energy just makes stuff too complicated
        double GB_energy; //Gay-Berne potential energy of the system that would be lost if this particle disappeared
        double bond_energy; //Energy associated with any bonds attached to this disk
        double wlc_energy; //Energy associated with any wlc bond attached to this disk
        double tot_energy; //All energy associated with this disk
        */


};
