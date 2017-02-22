#include "Disk.hpp"
#include "Box.hpp"
#include "mersenne.h"
#include "defs.h"

class Sim {
	public:
        Box box;
        Disk *disk;
        CRandomMersenne *RanGen;

        Sim();
        void MCSim();
        void openFiles();
        void read_input();
        void initialize_data();
        void initialize_system();
        void allocate_memory();
        void generate_chains();
        void calc_random_vector(Vector3d &b);
        void calc_random_normal_vector(Vector3d &b, Vector3d a);
        void calc_cross_vector(Vector3d &c, Vector3d a, Vector3d b);
        double calc_GB_total();
        void update_disk_energies();
        double calc_disk_energy(int);
        double calc_disk_interval_energy(int,int);
        double calc_disk_GB_energy(int);
        double calc_disk_GB_energy_without(int,int);
        double calc_bond_total();
        double calc_wlc_total();
        double calc_GB(Disk d1,Disk d2);
        double eps_func(double,double,double);
        double sig_func(double,double,double);
        double calc_bond(Disk d1, Disk d2);
        double calc_wlc(Disk d1, Disk d2);
        bool withinExclude(int,int);
        double calc_dist(Vector3d r1, Vector3d r2);
        void nearest_image_dist(Vector3d &r, Vector3d r1, Vector3d r2);
        void PBC_shift(Vector3d &r_shift, Vector3d r);
        void adjust_single_u_vector(int);
        void adjust_u_vectors(int);
        void proj_vector(Vector3d &a, Vector3d n);

        void MCMoveContin();

        void MCMove();
        void MC_displace();
        void MC_rotate();
        void MC_bend();
        void MC_twist();
        void MC_curl();
        void MC_reptate();
        void MC_translate();

        void shiftCOM();
        void printXYZ();
        void printPSF();
        void printPOV();
        void writeEnergy();
        void dumpData();

        void printCheckpoint();
        void readCheckpoint();
};
