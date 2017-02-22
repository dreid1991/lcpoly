#include "math_vector.h"
#include "Grid.h"
#include "defs.h"

using namespace LAMMPS_NS;

class Box {
	public:
       double x; //Dimensions of the periodic domain
       double y;
       double z;

       double T; //Temperature used for Monte Carlo simulation
       double E_tot, E_nb, E_bond, E_wlc; //Total energy, non-bond energy, harmonic bond energy, twistable worm-like chain energy
       int seed; //Random number seed
       int numDisks; //Number of disks on a single chain
       int numChains; //Number of chains total in the simulation
       int numTotal; //Total disks in the simulation

       double eps0,sig0,kappa,kappa_prime,chi,chi_prime,nu,mu; //Numbers used for calculating Gay-Berne interaction
       
       double kappa_b1, kappa_b2, kappa_t; //Parameters describing cost of bending and twisting for twistable worm like chain.

       double k_spr,r0; //Parameters to describe energy between harmonic bonds

       int excludeLength; //Distance along chain to exclude any non-bonded interactions

       double r_cut;

       int numCycles; //Number of Monte Carlo Cycles
       int cycle; //Current cycle
       int numAccepts;
       int numRejects;
       int vizInterval; //Interval after which visualization data is outputted
       int checkpoint; //0 if starting from scratch. 1 if loading a configuration
       int checkInterval; //MC cycle interval after which checkpoints are made

       double fracDisplace, fracRotate, fracBend, fracTwist, fracReptate, fracCurl, fracTranslate; //Relative occurence of the different Monte Carlo moves
       int maxReptate;
       int mode; // 0->particle, 1->contin
       
       //continuum stuff
       Grid<double> densities[2];
       Grid<Tensor> directorOrdering;


};