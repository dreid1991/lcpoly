#include "time.h"
#include <iostream>
#include <math.h>
#include <cstdlib>
#include <stdio.h>
#include "Sim.hpp"
#include "mersenne.h"
#define PI 3.14159265358979323846


enum MODE {PARTICLE, CONTINUUM};

using namespace LAMMPS_NS;
//using namespace std;

Sim::Sim() {
    read_input();
    initialize_data();
    initialize_system();
    shiftCOM();
    openFiles();
}

void Sim::MCSim() {
    box.numAccepts = 0;
    box.numRejects = 0;
    for (box.cycle=0;box.cycle<box.numCycles;box.cycle++) {
        if (box.cycle%box.vizInterval == 0) {
            shiftCOM();
            printXYZ();
            dumpData();
            writeEnergy();
            printf("Turn: %d\n",box.cycle);
        }
        if (box.mode==MODE::PARTICLE) {
            MCMove();
        } else {
            MCMoveContin();
        }
        if (box.cycle%box.checkInterval == 0) {
            printCheckpoint();
        }
    }
    writeEnergy();
    dumpData();
    printXYZ();
    printPOV();
    printf("Percent Acceptance: %f\n",double(box.numAccepts)/double(box.numAccepts+box.numRejects));
}

void Sim::openFiles() {
    FILE *dump;
    dump = fopen("config.xyz","w");
    fclose(dump);
    dump=fopen("Energy.txt","w");
    fclose(dump);
    dump = fopen("data.txt","w");
    fprintf(dump,"%d\n%d\n",box.numChains,box.numDisks);
    fclose(dump);
    printPSF();

}

void Sim::read_input() {

  FILE *input;
  char tt[2001];
  int ret;
  input=fopen("input", "r");

  if (input!=NULL)
    { 	/******		Reading Input File		********/
      ret=fscanf(input,"%d", &box.seed);		     fgets(tt,2000,input);
      ret=fscanf(input,"%d", &box.checkpoint);		 fgets(tt,2000,input);
      ret=fscanf(input,"%d", &box.numCycles);		 fgets(tt,2000,input);
      ret=fscanf(input,"%d", &box.vizInterval);		 fgets(tt,2000,input);
      ret=fscanf(input,"%d", &box.checkInterval);    fgets(tt,2000,input);
      ret=fscanf(input,"%lf", &box.fracDisplace);    fgets(tt,2000,input);
      ret=fscanf(input,"%lf", &box.fracRotate);		 fgets(tt,2000,input);
      ret=fscanf(input,"%lf", &box.fracBend);		 fgets(tt,2000,input);
      ret=fscanf(input,"%lf", &box.fracTwist);		 fgets(tt,2000,input);
      ret=fscanf(input,"%lf", &box.fracCurl);		 fgets(tt,2000,input);
      ret=fscanf(input,"%lf", &box.fracTranslate);   fgets(tt,2000,input);
      ret=fscanf(input,"%lf", &box.x);               fgets(tt,2000,input);
      ret=fscanf(input,"%lf", &box.y);               fgets(tt,2000,input);
      ret=fscanf(input,"%lf", &box.z);               fgets(tt,2000,input);
      ret=fscanf(input,"%d", &box.numDisks);         fgets(tt,2000,input);
      ret=fscanf(input,"%d", &box.numChains);        fgets(tt,2000,input);
      ret=fscanf(input,"%lf", &box.T);               fgets(tt,2000,input);
      ret=fscanf(input,"%lf", &box.eps0);            fgets(tt,2000,input);
      ret=fscanf(input,"%lf", &box.sig0);            fgets(tt,2000,input);
      ret=fscanf(input,"%lf", &box.kappa_b1);        fgets(tt,2000,input);
      ret=fscanf(input,"%lf", &box.kappa_b2);        fgets(tt,2000,input);
      ret=fscanf(input,"%lf", &box.kappa_t);         fgets(tt,2000,input);
      ret=fscanf(input,"%lf", &box.k_spr);           fgets(tt,2000,input);
      ret=fscanf(input,"%lf", &box.r0);              fgets(tt,2000,input);
      ret=fscanf(input,"%lf", &box.r_cut);           fgets(tt,2000,input);
      ret=fscanf(input,"%d",  &box.excludeLength);   fgets(tt,2000,input);
      ret=fscanf(input,"%d",  &box.mode);           fgets(tt,2000,input);
      fgets(tt,2000,input);fgets(tt,2000,input);

      fclose(input);
      /***************************************/
    }

  else
    { 
      fprintf(stdout,"\n Input File Not Found\n");
      exit(EXIT_FAILURE);
    }
    


}

void Sim::initialize_data() {
    box.nu = 1.0;
    box.mu = 2.0;
    box.kappa = 1.0/3;
    box.kappa_prime = 1.0/5;
    box.chi = (box.kappa*box.kappa-1)/(box.kappa*box.kappa+1);
    box.chi_prime = (sqrt(box.kappa_prime)-1)/(sqrt(box.kappa_prime)+1);
    box.numTotal = box.numDisks*box.numChains;
    if (box.seed < 0) {
        box.seed = time(NULL);
    }
    RanGen = new CRandomMersenne(box.seed);
    
}

void Sim::initialize_system() {
    double E_nb,E_bond,E_wlc;

    if (box.checkpoint == 1 && fopen("traj.bin","rb")!=NULL) {
        readCheckpoint();
    } else {
        box.checkpoint = 0;
        allocate_memory();
        generate_chains();
    }
    E_nb = calc_GB_total();
    E_bond = calc_bond_total();
    E_wlc = calc_wlc_total();
    box.E_tot = E_nb + E_bond + E_wlc;
}

void Sim::allocate_memory() {
    disk = new Disk[box.numTotal]; 
}

void Sim::generate_chains() {
    int i,j,k,m,m0;
    double x,y,z,R1,R2,R3;

    for (i=0;i<box.numChains;i++) {
        for (j=0;j<box.numDisks;j++) {
            m = i*box.numDisks+j;
            if (j==0) { //First disk in chain has random placement and random orientation 

                x = box.x * RanGen->Random();
                y = box.y * RanGen->Random();
                z = box.z * RanGen->Random();

                disk[m].rn[0] = x;
                disk[m].rn[1] = y;
                disk[m].rn[2] = z;

                PBC_shift(disk[m].r, disk[m].rn);

                calc_random_vector(disk[m].u); //Generate u completely randomly
                calc_random_normal_vector(disk[m].f, disk[m].u); //f is a vector normal to u 
                calc_cross_vector(disk[m].v, disk[m].u,disk[m].f); //u x f = v

                
            } else {//Other disks are placed in a perfectly straight line     
                for (k=0;k<3;k++) {
                    disk[m].u[k] = disk[m-1].u[k];
                    disk[m].f[k] = disk[m-1].f[k];
                    disk[m].v[k] = disk[m-1].v[k];

                    disk[m].rn[k] = disk[m-1].rn[k] + box.r0*disk[m-1].u[k]; 
                }
                
                PBC_shift(disk[m].r,disk[m].rn);
                
            }
        }
    }
}

void Sim::calc_random_vector(Vector3d &b) {
    double R1,R2,R3;
    do { //Generate random unit vector and make it u.
        R1 = (2*RanGen->Random()-1);
        R2 = (2*RanGen->Random()-1);
        R3 = R1*R1+R2*R2;
    } while (R3>=1);
               
    b[0] = 2*sqrtl(1.0-R3)*R1;
    b[1] = 2*sqrtl(1.0-R3)*R2;
    b[2] = 1-2.0*R3;


}

void Sim::calc_random_normal_vector(Vector3d &n, Vector3d a) { //Uses a rotation matrix to calculate a unit vector, n, normal to input vector a
    double theta, phi, A,B,C,D,R1,x0,y0;
    double ax,ay,az;
    ax = a[0]; ay = a[1]; az = a[2];
    theta = acos(az / sqrt(ax*ax+ay*ay+az*az));
    phi = atan2(ay,ax);
    A = cos(theta);
    B = sin(theta);
    C = cos(phi);
    D = sin(phi);
    R1 = 2*PI*RanGen->Random();
    x0 = cos(R1);
    y0 = sin(R1);
    n[0] = A*C*x0 - D*y0;
    n[1] = A*D*x0 + C*y0;
    n[2] = -1*B*x0;

    
}

void Sim::calc_cross_vector(Vector3d &c, Vector3d a, Vector3d b) {//Cross product
    
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];


}

/*********************************ENERGY CALCULATIONS**************************************/

double Sim::calc_GB_total() { // Does not double count the energy
    int i,j;
    double energy = 0.0;
    for (i=0;i<box.numTotal;i++) {
        for (j=i+1;j<box.numTotal;j++) {
            if (!withinExclude(i,j)) {
                energy += calc_GB(disk[i],disk[j]);
            }
        }
    }
    return energy;
}

double Sim::calc_disk_energy(int i) { //Computes all energy associated with a single disk (GB, bond, wlc) and updates its values
    int flag;
    double wlc = 0.0;
    double bond = 0.0;
    double GB;
    double tot;
    GB = calc_disk_GB_energy(i);

    if (i%box.numDisks==0) {
        flag = 0;
    } else if ((i+1)%box.numDisks==0) {
        flag = 2;
    } else {
        flag = 1;
    }

    if (flag!=2) {
        bond += calc_bond(disk[i],disk[i+1]);
        wlc += calc_wlc(disk[i],disk[i+1]);
    }
    if (flag!=0) {
        bond += calc_bond(disk[i-1],disk[i]);
        wlc += calc_wlc(disk[i-1],disk[i]);
    }

    tot = GB+bond+wlc;
    
    /* Maybe one day
    disk[i].GB_energy = GB;
    disk[i].bond_energy = bond;
    disk[i].wlc_energy = wlc;
    disk[i].tot_energy = tot;
    */

    return tot;

}

double Sim::calc_disk_interval_energy(int i, int j) { //Calculate all the energy (GB, bond, wlc) associated with the chain segment between disk i and j (i<j). Stuff gets complicated for multiple pair potentials
    int k;
    int flag1 = 0;
    int flag2 = 0;
    double wlc = 0.0;
    double bond = 0.0;
    double GB = 0.0;
    double tot;
    if (i%box.numDisks==0) { //First disk is start of chain
        flag1 = 1;
    } if ((j+1)%box.numDisks==0) { //Last disk is end of chain
        flag2 = 1;
    } 

    if (flag1 == 0) {
        bond += calc_bond(disk[i-1],disk[i]);
        wlc += calc_wlc(disk[i-1],disk[i]);
    }

    for (k=i;k<j;k++) {
        GB += calc_disk_GB_energy_without(k,i); //Prevents double counting of energies
        bond += calc_bond(disk[k],disk[k+1]);
        wlc += calc_wlc(disk[k],disk[k+1]);
    }
    GB += calc_disk_GB_energy_without(j,i);
    if (flag2==0) {
        bond += calc_bond(disk[j],disk[j+1]);
        wlc += calc_wlc(disk[j],disk[j+1]);
    }

    tot = GB+bond+wlc;
    return tot;

}

double Sim::calc_disk_GB_energy(int i) { //Calculates the Gay-Berne energy associated with a single disk
    int j;
    double energy = 0.0;
    for (j=0;j<box.numTotal;j++) {
        if(!withinExclude(i,j))
            energy += calc_GB(disk[i],disk[j]);
    }
    return energy;
}

double Sim::calc_disk_GB_energy_without(int i, int start_exclude) { //Calculates Gay-Berne energy associated with a single disk minus the disks whose energies have been tabulated -- start_exclude to i-1
    int j;
    double energy = 0.0;
    for (j=0;j<box.numTotal;j++) {
        if(!withinExclude(i,j) && (j < start_exclude || j > i)) {
            energy += calc_GB(disk[i],disk[j]);
        }
    }
    return energy;
}

double Sim::calc_bond_total() {
    int i,j,m;
    double energy = 0.0;
    for (i=0;i<box.numChains;i++) {
        for (j=0;j<box.numDisks-1;j++) {
            m = i*box.numDisks+j;
            energy += calc_bond(disk[m],disk[m+1]);

        }
    }
    return energy;
}

double Sim::calc_wlc_total() {
    int i,j,m;
    double energy = 0.0;
    for (i=0;i<box.numChains;i++) {
        for (j=0;j<box.numDisks-1;j++) {
            m = i*box.numDisks+j;
            energy += calc_wlc(disk[m],disk[m+1]);
        }
    }
    return energy;
}

double Sim::calc_GB(Disk d1, Disk d2) { //Gay-Berne energy is calculated between disk i and j
    int ii;
    Vector3d r;
    double f0,f1,f2,ener,r_dist;
    double sig,eps,a,a3,a6;
    nearest_image_dist(r, d1.r,d2.r);
    r_dist = sqrt(r.dot(r));

    if (r_dist > box.r_cut) {
        //printf("%f\n",ener);
        return 0.0;
    }

    f0 = d1.f.dot(d2.f);
    f1 = r.dot(d1.f);
    f2 = r.dot(d2.f);
    eps = eps_func(f0,f1,f2);
    sig = sig_func(f0,f1,f2);
    a = box.sig0 / (r_dist - sig + box.sig0);
    a3 = a*a*a;
    a6 = a3*a3;

    ener = 4*eps*(a6*(a6-1));
    //printf("%f\n",ener);
    return ener;
}

double Sim::eps_func(double f0, double f1, double f2) {
    double v1,v2,v3;
    double out;

    v1 = 1.0/sqrt(1-box.chi*box.chi*f0*f0);
    v2 = (f1+f2)*(f1+f2)/(1+box.chi_prime*f0) + (f1-f2)*(f1-f2)/(1-box.chi_prime*f0);
    v3 = 1-box.chi_prime*v2/2;

    out = box.eps0 * v1*v3*v3;
    return out;

}

double Sim::sig_func(double f0, double f1, double f2) {
    double v1,out;

    v1 = (f1+f2)*(f1+f2)/(1+box.chi*f0) + (f1-f2)*(f1-f2)/(1-box.chi*f0);
    out = box.sig0 / sqrt(1 - box.chi*v1/2);
    return out;    

}

double Sim::calc_bond(Disk d1, Disk d2) { //Energy of bond is calculated between disk i and i+1
    double r;
    double e;
    Vector3d v1 = {d1.rn[0],d1.rn[1],d1.rn[2]};
    Vector3d v2 = {d2.rn[0],d2.rn[1],d2.rn[2]};
    r = calc_dist(v1,v2);
    e = 0.5*box.k_spr*(r-box.r0)*(r-box.r0);
    return e;
}

double Sim::calc_wlc(Disk d1, Disk d2) { //Twisting and bending energy is calculated between disk i and i+1
    int j;
    double w1sqr,w2sqr,w3sqr,theta,value;
    double du,df,dv;
    double e;
    du = calc_dist(d1.u,d2.u);
    df = calc_dist(d1.f,d2.f);
    dv = calc_dist(d1.v,d2.v);
    du *= du;
    df *= df;
    dv *= dv;
    w1sqr = 0.5*(du+dv-df);
    w2sqr = 0.5*(du+df-dv);
    w3sqr = 0.5*(df+dv-du);

    //Calculation of dihedral angle between disk 1 and disk 2

    Vector3d n, v1, v2; //Normal vector n defined by the vector between disks, v1 and v2 will be disk 1 and 2's normal (f) vector projected onto plane normal to n from which dihedral angle can be calculated
    for (j=0;j<3;j++) {
        n[j] = d2.rn[j] - d1.rn[j];
        v1[j] = d1.f[j];
        v2[j] = d2.f[j];
    }
    n.normalize();
    proj_vector(v1,n);
    proj_vector(v2,n);
    value = v1.dot(v2);
    if (value > 1.0) {
        value = 1.0;
    }
    if (value < -1.0) {
        value = -1.0;
    }
    theta = acos(value);

    //e = box.kappa_b1*w1sqr + box.kappa_b2*w2sqr + box.kappa_t*w3sqr; //For full twistable wormlike chain
    e = 0.5*box.kappa_b1*w1sqr + 0.5*box.kappa_b2*w2sqr + box.kappa_t*(1 - cos(2*theta));
    
    return e;
    
}


/*************************************************************************************************/

bool Sim::withinExclude(int i, int j) { //Returns true if i and j are on the same chain and within the exclude distance
    int chain1,chain2;
    chain1 = int(floor(double(i)/double(box.numDisks)));
    chain2 = int(floor(double(j)/double(box.numDisks)));
    if (chain1==chain2 && abs(i-j) <= box.excludeLength) {
        return true;
    } else {
        return false;
    }
}

double Sim::calc_dist(Vector3d r1, Vector3d r2) {

    double x,y,z,r;
    x = r1[0] - r2[0];
    y = r1[1] - r2[1];
    z = r1[2] - r2[2];

    r = sqrt(x*x+y*y+z*z);
    return r;

}

void Sim::nearest_image_dist(Vector3d &r, Vector3d r1, Vector3d r2) { //r1 - r2 using the nearest image convention
    double dx,dy,dz;
    double xhalf = box.x/2;
    double yhalf = box.y/2;
    double zhalf = box.z/2;

    dx = r1[0] - r2[0];
    dy = r1[1] - r2[1];
    dz = r1[2] - r2[2];

    dx = dx > xhalf ? dx - xhalf : dx; 
    dx = dx < -xhalf ? dx + xhalf : dx; 
    dy = dy > yhalf ? dy - yhalf : dy; 
    dy = dy < -yhalf ? dy + yhalf : dy; 
    dz = dz > zhalf ? dz - zhalf : dz; 
    dz = dz < -zhalf ? dz + zhalf : dz; 

    r[0] = dx; r[1] = dy; r[2] = dz;

}

void Sim::PBC_shift(Vector3d &r_shift, Vector3d r) {
    r_shift[0] = r[0] - box.x*floor(r[0]/box.x);
    r_shift[1] = r[1] - box.x*floor(r[1]/box.y);
    r_shift[2] = r[2] - box.x*floor(r[2]/box.z);
}

void Sim::adjust_single_u_vector(int i) {
        if (i%box.numDisks==0) { //Change the u vector so that it remains tangent to the direction of the chain
            disk[i].u[0] = disk[i+1].rn[0] - disk[i].rn[0];
            disk[i].u[1] = disk[i+1].rn[1] - disk[i].rn[1];
            disk[i].u[2] = disk[i+1].rn[2] - disk[i].rn[2];

            disk[i].u.normalize();
            proj_vector(disk[i].f,disk[i].u);
            disk[i].v = disk[i].u.cross(disk[i].f);
        } else if ((i+1)%box.numDisks==0) {
            disk[i].u[0] = disk[i].rn[0] - disk[i-1].rn[0];
            disk[i].u[1] = disk[i].rn[1] - disk[i-1].rn[1];
            disk[i].u[2] = disk[i].rn[2] - disk[i-1].rn[2];

            disk[i].u.normalize();
            proj_vector(disk[i].f,disk[i].u);
            disk[i].v = disk[i].u.cross(disk[i].f);

        } else {
            disk[i].u[0] = (disk[i+1].rn[0]-disk[i-1].rn[0])/2;
            disk[i].u[1] = (disk[i+1].rn[1]-disk[i-1].rn[1])/2;
            disk[i].u[2] = (disk[i+1].rn[2]-disk[i-1].rn[2])/2;
            disk[i].u.normalize();
            proj_vector(disk[i].f,disk[i].u);
            disk[i].v = disk[i].u.cross(disk[i].f);

        }

}

void Sim::adjust_u_vectors(int i) { //Some of the moves require a correction to the u vector to make sure that it is tangent to the chain. This function corrects it to make sure that it works 
                                    //WARNING: Will not work for chains with only two disks
        if (i%box.numDisks==0) { 
            disk[i].u[0] = disk[i+1].rn[0] - disk[i].rn[0];
            disk[i].u[1] = disk[i+1].rn[1] - disk[i].rn[1];
            disk[i].u[2] = disk[i+1].rn[2] - disk[i].rn[2];

            disk[i].u.normalize();
            proj_vector(disk[i].f,disk[i].u);
            disk[i].v = disk[i].u.cross(disk[i].f);

            disk[i+1].u[0] = (disk[i+2].rn[0] - disk[i].rn[0])/2;
            disk[i+1].u[1] = (disk[i+2].rn[1] - disk[i].rn[1])/2;
            disk[i+1].u[2] = (disk[i+2].rn[2] - disk[i].rn[2])/2;

            disk[i+1].u.normalize();
            proj_vector(disk[i+1].f,disk[i+1].u);
            disk[i+1].v = disk[i+1].u.cross(disk[i+1].f);

        } else if ((i+1)%box.numDisks==0) {
            disk[i].u[0] = disk[i].rn[0] - disk[i-1].rn[0];
            disk[i].u[1] = disk[i].rn[1] - disk[i-1].rn[1];
            disk[i].u[2] = disk[i].rn[2] - disk[i-1].rn[2];

            disk[i].u.normalize();
            proj_vector(disk[i].f,disk[i].u);
            disk[i].v = disk[i].u.cross(disk[i].f);


            disk[i-1].u[0] = (disk[i].rn[0] - disk[i-2].rn[0])/2;
            disk[i-1].u[1] = (disk[i].rn[1] - disk[i-2].rn[1])/2;
            disk[i-1].u[2] = (disk[i].rn[2] - disk[i-2].rn[2])/2;

            disk[i-1].u.normalize();
            proj_vector(disk[i-1].f,disk[i-1].u);
            disk[i-1].v = disk[i-1].u.cross(disk[i-1].f);

        } else {
            disk[i+1].u[0] = (disk[i+2].rn[0] - disk[i].rn[0])/2;
            disk[i+1].u[1] = (disk[i+2].rn[1] - disk[i].rn[1])/2;
            disk[i+1].u[2] = (disk[i+2].rn[2] - disk[i].rn[2])/2;

            disk[i+1].u.normalize();
            proj_vector(disk[i+1].f,disk[i+1].u);
            disk[i+1].v = disk[i+1].u.cross(disk[i+1].f);


            disk[i-1].u[0] = (disk[i].rn[0] - disk[i-2].rn[0])/2;
            disk[i-1].u[1] = (disk[i].rn[1] - disk[i-2].rn[1])/2;
            disk[i-1].u[2] = (disk[i].rn[2] - disk[i-2].rn[2])/2;

            disk[i-1].u.normalize();
            proj_vector(disk[i-1].f,disk[i-1].u);
            disk[i-1].v = disk[i-1].u.cross(disk[i-1].f);


        }

}

void Sim::proj_vector(Vector3d &a, Vector3d n) {//Project vector a onto the plane described by unit normal vector n, then normalize the new a
    double dot;
    int i;
    Vector3d n2;
    for (i=0;i<3;i++) {
        n2[i] = n[i];
    }
    dot = a.dot(n2);
    for (i=0;i<3;i++) {
        a[i] = a[i] - dot*n2[i];
    }
    a.normalize();
}




void Sim::MCMoveContin() {
    
}

/***************************************************MONTE CARLO MOVES******************************************/

void Sim::MCMove() {
    double displace = box.fracDisplace;
    double rotate = box.fracDisplace+box.fracRotate;
    double bend = box.fracDisplace+box.fracRotate+box.fracBend;
    double twist = box.fracDisplace+box.fracRotate+box.fracBend+box.fracTwist;
    double curl = box.fracDisplace+box.fracRotate+box.fracBend+box.fracTwist+box.fracCurl;
    double translate = box.fracDisplace+box.fracRotate+box.fracBend+box.fracTwist+box.fracCurl+box.fracTranslate;
    double rand = RanGen->Random();
    
    if (rand < displace) {
      //  MC_displace();
    } else if (rand < rotate) {
      //  MC_rotate();
    } else if (rand < bend) {
      //  MC_bend();
    } else if (rand < twist) {
        MC_twist();
    } else if (rand < curl) {
        //MC_curl();
    } else if (rand < translate) {
        //MC_translate();
    }
}

void Sim::MC_displace() { //Displace a disk slightly without changing its orientation
    int i,j,k,flag,first,last;
    double rx,ry,rz;
    double delta = box.r0/2; //Maximum displacement distance in any direction
    double pre_energy, post_energy, dE;
    Disk *save = new Disk[3];
    
    for (int move=0;move<box.numTotal;move++) {
        do {
            i = double(box.numTotal)*RanGen->Random();
        } while(i >= box.numTotal);
        
        rx = (2.0*RanGen->Random()-1.0)*delta;
        ry = (2.0*RanGen->Random()-1.0)*delta;
        rz = (2.0*RanGen->Random()-1.0)*delta;

        if (i%box.numDisks==0) {
            first = i; last = i+1;
        } else if ((i+1)%box.numDisks==0) {
            first = i-1; last = i;
        } else {
            first = i-1; last = i+1;
        }
        k=0;
        for (j=first; j<=last;j++) {
            save[k] = disk[j];
            k++;
        }
        
        pre_energy = calc_disk_interval_energy(first,last);
        
        disk[i].rn[0] += rx;
        disk[i].rn[1] += ry;
        disk[i].rn[2] += rz;
        
        PBC_shift(disk[i].r, disk[i].rn);
      
        adjust_u_vectors(i);

        post_energy = calc_disk_interval_energy(first,last);

        dE = post_energy - pre_energy;
        if (RanGen->Random() < exp(-dE/box.T)) { //Accept changes
            box.E_tot += dE;
            box.numAccepts++;
        } else { //Reverse changes
            k = 0;
            for (j=first;j<=last;j++) {
                disk[j] = save[k];
                k++;
            }
            box.numRejects++;
        }
         
        
    }

    delete [] save;
}

void Sim::MC_rotate() { //Randomly rotate a disk about some axis without changing its center of mass
    int i,j,flag;
    double max_theta = PI/2; //Don't want more than a 90 degree rotation around any axis
    double theta;
    double pre_energy, post_energy, dE;
    Disk save;
    Vector3d axis;

    for (int move=0;move<box.numTotal;move++) {
        do {
            i = double(box.numTotal)*RanGen->Random();
        } while(i >= box.numTotal);
        
        theta = max_theta*(2*RanGen->Random()-1); //Random amount to rotate by
       
        pre_energy = calc_disk_energy(i);
        save = disk[i];

        AngleAxisd aaU(theta, save.u);
        Matrix3d rot; 
        rot = aaU;
        disk[i].f = rot * save.f;
        disk[i].v = rot * save.v;


        post_energy = calc_disk_energy(i);

        dE = post_energy - pre_energy;

        if (RanGen->Random() < exp(-dE/box.T)) { //Accept changes
            box.E_tot += dE; 
            box.numAccepts++;
        } else { //Reverse changes
            disk[i] = save;
            box.numRejects++;
        }
           

    }
}

void Sim::MC_bend() { //Beginning at some point along the backbone, bend the chain uniformly (keep orientations of bending section the same relative to one another)
    int rand_chain,i,j,k,flag,dir,core,disk_length,ref,first,last;
    double max_theta = PI/3; //Don't want more than a 60 degree rotation around any axis
    double theta;
    double pre_energy, post_energy, dE;
    Vector3d axis; // Axis of rotation
    Vector3d dist1; // Contains distance vector for each disk from center of mass for rotation
    Vector3d dist2; // Contains new distance vector for each disk from center of mass for rotation
    Vector4d quat;
    Disk *chain = new Disk[box.numDisks];

    for (int move=0;move<box.numChains;move++) {
        
        do {
            rand_chain = double(box.numChains)*RanGen->Random();
        } while(rand_chain >= box.numChains);
        
        do {
            core = double(box.numDisks)*RanGen->Random();
        } while(core >= box.numDisks);
        
        ref = rand_chain*box.numDisks+core; 

        if (2*RanGen->Random()-1 > 0) {
            dir = 1; //
            first = ref;
            last = box.numDisks*rand_chain + box.numDisks-1; //last disk in current chain
        } else {
            first = box.numDisks*rand_chain; //first disk in current chain
            last = ref;
            dir = -1;
        }
        pre_energy = calc_disk_interval_energy(first,last);
        j = 0;
        for (i=first;i<=last;i++) {
            chain[j] = disk[i];
            j++;
        }
        

        /*
            AngleAxisd aa(theta, disk[i].f);
            Matrix3d rot;
            rot = aa;
            Vector3d v = rot * disk[i].u;
            */


        theta = max_theta*(2*RanGen->Random()-1); //Random amount to rotate by
        calc_random_vector(axis); //Random unit vector for rotation
        //quat[0] = cos(theta/2);
        for (j=0;j<3;j++) {
          //  quat[j+1] = axis[j]*sin(theta/2);
        }
       
        AngleAxisd aa(theta, axis);
        Matrix3d rot;
        rot = aa;
        j = 0;
        for (i=first;i<=last;i++) {
            for (k=0;k<3;k++) {
                dist1[k] = disk[i].rn[k] - disk[ref].rn[k];
            }
            //quat_vec_rot(dist2,dist1,quat);
            dist2 = rot * dist1;
            for (k=0;k<3;k++) {
                disk[i].rn[k] = disk[ref].rn[k] + dist2[k];
            }
            PBC_shift(disk[i].r, disk[i].rn);
            disk[i].u = rot * chain[j].u;
            disk[i].f = rot * chain[j].f;
            disk[i].v = rot * chain[j].v;
            //quat_vec_rot(disk[i].u,chain[j].u,quat);
            //quat_vec_rot(disk[i].f,chain[j].f,quat);
            //quat_vec_rot(disk[i].v,chain[j].v,quat);
            j++;

        }

        adjust_single_u_vector(ref);

        post_energy = calc_disk_interval_energy(first,last);

        dE = post_energy - pre_energy;

        if (RanGen->Random() < exp(-dE/box.T)) { //Accept changes
            box.E_tot += dE; 
            box.numAccepts++;
        } else { //Reverse changes
            j=0;
            for (i=first;i<=last;i++) {
                disk[i] = chain[j];
                j++;
            }
            box.numRejects++;
        }
    }
    
    delete [] chain;
}

void Sim::MC_twist() { //Beginning at some point along the backbone, uniformly twist the chain
    int rand_chain,i,j,k,flag,dir,core,disk_length,ref,first,last,increment;
    double max_dtheta = PI/12; //Keeps incremental twists below 15 degrees
    double theta, dtheta;
    double pre_energy, post_energy, dE;
    Vector4d quat;
    Disk *chain = new Disk[box.numDisks];

    for (int move=0;move<box.numChains;move++) {
        
        do {
            rand_chain = double(box.numChains)*RanGen->Random();
        } while(rand_chain >= box.numChains);
        
        do {
            core = double(box.numDisks)*RanGen->Random();
        } while(core >= box.numDisks);
        
        ref = rand_chain*box.numDisks+core; 

        if (2*RanGen->Random()-1 > 0) {
            dir = 1; //
            first = ref;
            last = box.numDisks*rand_chain + box.numDisks-1; //last disk in current chain
        } else {
            first = box.numDisks*rand_chain;
            last = ref;
            dir = -1;
        }

        pre_energy = calc_disk_interval_energy(first,last);
        j = 0;
        for (i=first;i<=last;i++) {
            chain[j] = disk[i];
            j++;
        }
        
        dtheta = max_dtheta*RanGen->Random(); //Random amount to rotate by
       
        j = 0;
        if (dir == 1) increment = 1;
        if (dir == -1) increment = last-first+1;
        for (i=first;i<=last;i++) {
            
            theta = (increment)*dtheta;
            AngleAxisd aa(theta, disk[i].u);
            Matrix3d rot;
            rot = aa;
            //quat[0] = cos(theta/2);
            for (k=0;k<3;k++) {
              //  quat[k+1] = disk[i].u[k]*sin(theta/2);
            }
            disk[i].f = rot * chain[j].f;
            disk[i].v = rot * chain[j].v;
            //quat_vec_rot(disk[i].f,chain[j].f,quat);
            //quat_vec_rot(disk[i].v,chain[j].v,quat);
            j++;
            increment += dir;

        }

        post_energy = calc_disk_interval_energy(first,last);

        dE = post_energy - pre_energy;

        if (RanGen->Random() < exp(-dE/box.T)) { //Accept changes
            box.E_tot += dE; 
            box.numAccepts++;
        } else { //Reverse changes
            j=0;
            for (i=first;i<=last;i++) {
                disk[i] = chain[j];
                j++;
            }
            box.numRejects++;
        }
    }

    delete [] chain;

}

void Sim::MC_curl() { //Imparts curvature to the chain by "curling" it in a randomly chosen direction
    int rand_chain,i,j,k,flag,dir,core,disk_length,ref,first,last,tot_move;
    double max_theta = PI/50.0; //Low theta values because it can result in a lot of curvature  
    double theta,tot_theta;
    double pre_energy, post_energy, dE;
    Vector3d axis; // Axis of rotation
    Vector3d dist1; // Contains distance Vector3d for each disk from center of mass for rotation
    Vector3d dist2; // Contains new distance vector for each disk from center of mass for rotation
    Quaterniond quat;
    Disk *chain = new Disk[box.numDisks];

    for (int move=0;move<box.numChains;move++) {
        
        do {
            rand_chain = double(box.numChains)*RanGen->Random();
        } while(rand_chain >= box.numChains);
        
        do {
            core = double(box.numDisks)*RanGen->Random();
        } while(core >= box.numDisks);
        
        ref = rand_chain*box.numDisks+core; 

        if (2*RanGen->Random()-1 > 0) {
            dir = 1; //
            first = ref;
            last = box.numDisks*rand_chain + box.numDisks-1; //last disk in current chain
            core = 0;
        } else {
            first = box.numDisks*rand_chain; //first disk in current chain
            last = ref;
            dir = -1;
            core = last-first;
        }
        tot_move = last-first+1;
        pre_energy = calc_disk_interval_energy(first,last);
        j = 0;
        for (i=first;i<=last;i++) {
            chain[j] = disk[i];
            j++;
        }
        
        theta = max_theta*(2*RanGen->Random()-1); //Random amount to rotate by
        calc_random_normal_vector(axis,disk[ref].u); //Random unit vector for rotation
       
        for (i=1;i<tot_move;i++) {
            j=i*dir;
            tot_theta = double(i)*theta;
            //quat[0] = cos(tot_theta/2);
            for (k=0;k<3;k++) {
                //quat[k+1] = axis[k]*sin(tot_theta/2);
                dist1[k] = chain[core+j].rn[k] - chain[core+j-dir].rn[k];
            }
            //quat_vec_rot(dist2,dist1,quat);
            for (k=0;k<3;k++) {
                disk[ref+j].rn[k] = disk[ref+j-dir].rn[k] + dist2[k];
            }
            PBC_shift(disk[ref+j].r, disk[ref+j].rn);
            //quat_vec_rot(disk[ref+j].u,chain[core+j].u,quat);
            //quat_vec_rot(disk[ref+j].f,chain[core+j].f,quat);
            //quat_vec_rot(disk[ref+j].v,chain[core+j].v,quat);

        }
        
        for (i=first;i<=last;i++) {
            adjust_single_u_vector(i);
        }

        post_energy = calc_disk_interval_energy(first,last);

        dE = post_energy - pre_energy;

        if (RanGen->Random() < exp(-dE/box.T)) { //Accept changes
            box.E_tot += dE; 
            box.numAccepts++;
        } else { //Reverse changes
            j=0;
            for (i=first;i<=last;i++) {
                disk[i] = chain[j];
                j++;
            }
            box.numRejects++;
        }
    }
    
    delete [] chain;
  
}

void Sim::MC_reptate() { //Snake-like slithering move

}

void Sim::MC_translate() { //Randomly translate the entire polymer -- This is made as a global move
    int rand_chain,i,j,k,flag,dir,core,disk_length,ref,first,last,tot_move;
    double x,y,z,xCenter,yCenter,zCenter;   
    double pre_energy, post_energy, dE;
    Vector3d dist1; // Contains distance Vector3d for each disk from center of mass for rotation
    Vector3d dist2; // Contains new distance Vector3d for each disk from center of mass for rotation
    Disk *chain = new Disk[box.numDisks];

    for (int move=0;move<box.numChains;move++) {
        
        do {
            rand_chain = double(box.numChains)*RanGen->Random();
        } while(rand_chain >= box.numChains);
       
        first = rand_chain*box.numDisks;
        last = (rand_chain+1)*box.numDisks-1;

        pre_energy = calc_disk_interval_energy(first,last);
        for (i=0;i<box.numDisks;i++) {
            chain[i] = disk[first+i];
        }
       
        xCenter = 0.0;
        yCenter = 0.0;
        zCenter = 0.0;
        for (j=first;j<=last;j++) {
            xCenter += disk[j].rn[0];
            yCenter += disk[j].rn[1];
            zCenter += disk[j].rn[2];
        }
        xCenter /= double(box.numDisks); //Calculation of center of mass of polymer
        yCenter /= double(box.numDisks);
        zCenter /= double(box.numDisks);
        x = box.x*RanGen->Random(); //Calculation for new center of mass
        y = box.y*RanGen->Random();
        z = box.z*RanGen->Random(); 
        x -= xCenter;               //Change in center of mass
        y -= yCenter;
        z -= zCenter;
        for (j=first;j<=last;j++) {
            disk[j].rn[0] += x;
            disk[j].rn[1] += y;
            disk[j].rn[2] += z;
            PBC_shift(disk[j].r,disk[j].rn);
        }

        post_energy = calc_disk_interval_energy(first,last);

        dE = post_energy - pre_energy;

        if (RanGen->Random() < exp(-dE/box.T)) { //Accept changes
            box.E_tot += dE; 
            box.numAccepts++;
        } else { //Reverse changes
            for (i=0;i<box.numDisks;i++) {
                disk[first+i] = chain[i];
            }
            box.numRejects++;
        }
    }
    
    delete [] chain;
  

}
/**************************************************************************************************************/


/**********************************DATA OUTPUT AND VISUALIZATION***********************************************/
void Sim::shiftCOM() { //Shifts the center of mass of all the nonperiodic coordinates back into the simulation box 
    int i,j,m;
    double xCenter,yCenter,zCenter,dx,dy,dz;
    for (i=0;i<box.numChains;i++) {
        xCenter = 0.0;
        yCenter = 0.0;
        zCenter = 0.0;
        for (j=0;j<box.numDisks;j++) {
            m = i*box.numDisks+j;
            xCenter += disk[m].rn[0];
            yCenter += disk[m].rn[1];
            zCenter += disk[m].rn[2];
        }
        xCenter /= double(box.numDisks);
        yCenter /= double(box.numDisks);
        zCenter /= double(box.numDisks);
        dx = box.x*floor(xCenter/box.x);
        dy = box.y*floor(yCenter/box.y);
        dz = box.z*floor(zCenter/box.z);
        for (j=0;j<box.numDisks;j++) {
            m = i*box.numDisks+j;
            disk[m].rn[0] -= dx;
            disk[m].rn[1] -= dy;
            disk[m].rn[2] -= dz;
        }
    }
}

void Sim::printXYZ() {
    int i,j;
    FILE *dump;
    Vector4d q;
    double dtheta = 2*PI/12;
    double theta;
    dump = fopen("config.xyz","a");
    fprintf(dump,"%d\nHELLO\n",13*box.numTotal);
    for (i=0;i<box.numTotal;i++) {
        fprintf(dump,"0\t%f\t%f\t%f\n",disk[i].rn[0],disk[i].rn[1],disk[i].rn[2]);
    }
    for (i=0;i<box.numTotal;i++) {
       for (j=0;j<12;j++) {
            theta = dtheta*j;
            AngleAxisd aa(theta, disk[i].f);
            Matrix3d rot;
            rot = aa;
            Vector3d v = rot * disk[i].u;

            fprintf(dump,"0\t%f\t%f\t%f\n",disk[i].rn[0]+v[0]*box.r0/3,disk[i].rn[1]+v[1]*box.r0/3,disk[i].rn[2]+v[2]*box.r0/3);
       }
    }

    fclose(dump);
}

void Sim::printPSF() {
      FILE *sout;
      int counter,type,pair,chain_switch,bead_type;
      char t[2]="O";
      sout=fopen("polymer.psf","w");
      fprintf(sout,"*\n*\n*\n*\n*\n\n");
      fprintf(sout,"%7d !NATOMS\n",13*box.numTotal);
      for(counter=1;counter<=13*box.numTotal;counter++){
        if(counter<=box.numTotal){

          {t[0]='O'; bead_type=1;}

          fprintf(sout, "%8d ", counter);
          fprintf(sout, "POLY ");
          fprintf(sout, "%-4d ",1);
          fprintf(sout, "%s  ", "POL");
          fprintf(sout, "%-5s ",t);
          fprintf(sout, "%3d  ", bead_type);
          fprintf(sout, "%13.6e   ",0.0);
          fprintf(sout, "%7.3lf           0\n", 1.0);
        }

        else{
          t[0]='S'; bead_type=2;
          fprintf(sout, "%8d ", counter);
          fprintf(sout, "SOLV ");
          fprintf(sout, "%-4d ",1);
          fprintf(sout, "%s  ", "SOL");
          fprintf(sout, "%-5s ",t);
          fprintf(sout, "%3d  ", bead_type);
          fprintf(sout, "%13.6e   ",0.0);
          fprintf(sout, "%7.3lf           0\n", 1.0);

        }        
        //fprintf(sout,"%8d\tU\t%4d\tPOL\t%1s\t%1s\t0.0000000\t1.0000\t0\n  ",counter,counter-1,t,t);
      }


      fprintf(sout,"\n%8d !NBOND: bond\n",(box.numDisks-1)*box.numChains);




      pair=1;
      counter=1;
      for(counter=1;counter<=box.numTotal;counter++){
        chain_switch=counter%box.numDisks;
        //printf("**chain_switch:\t%d\n",chain_switch);
        if(chain_switch==0)
          {counter++;
        if(counter>box.numTotal)
          break;
          }

        if(pair==5)
          {pair=1;fprintf(sout,"\n");}
        //printf("**counter:\t%d\n",counter);
        fprintf(sout,"%8d%8d",counter,counter+1);
        pair++;
      }

        fprintf(sout,"\n%8d !NTHETA: angles\n",0);
      fprintf(sout,"\n%8d !NPHI: dihedrals\n",0);
      fprintf(sout,"\n%8d !NIMPR\n",0);
      fprintf(sout,"\n%8d !HDON\n",0);
      fprintf(sout,"\n%8d !HACC\n",0);
      fprintf(sout,"\n%8d !NNB\n",0);


      fclose(sout);

}

void Sim::printPOV() {
    FILE *pov = fopen("disks.xyz","w");
    fprintf(pov, "ITEM: TURN\n%d\nITEM: NUMBER OF ATOMS\n%d\nITEM: BOX BOUNDS\n-10 10\n-10 10\n-10 10\nITEM: ATOMS\n",box.cycle,box.numTotal);

    double scalex = 10.0/(box.x/2);
    double scaley = 10.0/(box.y/2);
    double scalez = 10.0/(box.z/2);

    for (int i=0;i<box.numTotal;i++) {
        fprintf(pov,"%d 1 %f %f %f %f %f %f\n", i, (disk[i].r[0]-box.x/2)*scalex, (disk[i].r[1]-box.y/2)*scaley, (disk[i].r[2]-box.z/2)*scalez, disk[i].f[0], disk[i].f[1], disk[i].f[2]); 
    }
    fclose(pov);
}

void Sim::writeEnergy() {
    FILE *en = fopen("Energy.txt","a");
    fprintf(en,"%d\t%f\n",box.cycle,box.E_tot);
    fclose(en);
}

void Sim::dumpData() {
    FILE *out = fopen("data.txt","a");
    int i;
    fprintf(out,"%d\n",box.cycle);
    for (i=0;i<box.numTotal;i++) {
        fprintf(out,"%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",disk[i].rn[0],disk[i].rn[1],disk[i].rn[2],disk[i].u[0],disk[i].u[1],disk[i].u[2],disk[i].f[0],disk[i].f[1],disk[i].f[2],disk[i].v[0],disk[i].v[1],disk[i].v[2]);
    }
    fclose(out);
}
/**************************************************************************************************************/

/************************************TRAJECTORY CHECKPOINTING FUNCTIONS****************************************/

void Sim::printCheckpoint() {
    FILE *out = fopen("traj.bin","wb");
    int i;
    fwrite(&box, sizeof(Box), 1 , out);
    for (i=0;i<box.numTotal;i++) {
        fwrite(&disk[i], sizeof(Disk), 1 , out);
    }
    fclose(out);
}

void Sim::readCheckpoint() {
    FILE *in = fopen("traj.bin","rb");
    int i;
    fread(&box, sizeof(Box),1,in);
    allocate_memory();
    for (i=0;i<box.numTotal;i++) {
        fread(&disk[i],sizeof(Disk),1,in);
    }
    fclose(in);
}

/**************************************************************************************************************/
