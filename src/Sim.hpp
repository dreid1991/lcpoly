#include "Disk.hpp"
#include "Box.hpp"
#include "mersenne.h"
#include "defs.h"
#include <vector>
#include <utility>  //pair
#include "Polymer.h"

//typedef std::vector<std::pair<int, Disk> > StoredState; //vector of idx and disk pairs, for resetting after a rejected move
class Sim {
	public:
        Box box;
        std::vector<Disk> disk;
        std::vector<Polymer> poly;
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

        double calc_eng_contin_total();

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



        void set_bulk_density();
        void set_density();
        void add_density(Polymer &);
        void sub_density(Polymer &);
        void modify_density(Polymer &, int);

        void set_align_tensor_u();
        void add_align_tensor_u(Polymer &);
        void sub_align_tensor_u(Polymer &);
        void modify_align_tensor_u(Polymer &, int);
        
        void add_all_grids(Polymer &p);
        void sub_all_grids(Polymer &p);

        double calc_comp_total();
        double calc_chi_total();
        double calc_align_u_total();

        void MCMoveContin();
        void MC_contin_displace(); 
        void MC_contin_rotate(); 
        void MC_contin_translate();

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
