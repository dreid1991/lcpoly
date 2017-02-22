#include "math_vector.h"

using namespace LAMMPS_NS;

class Disk {
	public:
        vector r; //Contains x,y,z coordinates
        vector rn; //Contains x,y,z coordinates without periodic boundary conditions for bond calculation purposes
        vector u; //Vector tangent to the polymer 
        vector f; //Vector normal to the disk (used for Gay-Berne calculations)
        vector v; //Vector normal to u and f

        /* Using this to store energy just makes stuff too complicated
        double GB_energy; //Gay-Berne potential energy of the system that would be lost if this particle disappeared
        double bond_energy; //Energy associated with any bonds attached to this disk
        double wlc_energy; //Energy associated with any wlc bond attached to this disk
        double tot_energy; //All energy associated with this disk
        */


};
