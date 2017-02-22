/* This program is designed to calculate an energy histogram using dihedral energies and the tight binding model. It relies on Eigen for linear algebra calculations */
#include "time.h"
#include <iostream>
#include <math.h>
#include <cstdlib>
#include <stdio.h>
#define PI 3.14159265358979323846
#include "Eigen/Dense"
#include "math_vector.h"
#include "Disk.hpp"

using namespace std;
using namespace Eigen;
using namespace LAMMPS_NS;

int numChains,numDisks,numTotal,turn,numLines,numSnapshots; 
Disk *disk;
FILE *dataFile;
MatrixXd Hamiltonian;
MatrixXd Eigenvalues;

void init_data() {
    ifstream aFile("data.txt");
    numLines = 0;
    string line;
    while(getline(aFile,line)) 
        numLines++
    
    char tt[2001];
    dataFile = fopen("data.txt","r");
    fscanf(dataFile, "%d%d",&numChains,&numDisks); fgets(tt,2000,dataFile);
    numTotal = numDisks*numChains;
    disk = new Disk[numTotal];

    numSnapshots = (numLines - 2) /(numTotal+1); //Should be an exact integer
}

void init_matrix() {
    Hamiltonian = MatrixXd::Zero(numDisks,numDisks);
    Eigenvalues = MatrixXd::Zero(numDisks,numDisks);
}

void kill_everything() {
    fclose(dataFile);
    delete [] disk;
}

void diagonalize() {
    SelfAdjointEigenSolver<MatrixXd> solver(Hamiltonian);
    Eigenvalues = solver.eigenvalues().asDiagonal();
}

void readSnapshot() {
    int i;
    fscanf(dataFile,"%d",&turn);
    for (i=0;i<numTotal;i++) {
        fscanf(datafile,"%f%f%f%f%f%f%f%f%f%f%f%f",&disk[i].rn[0],&disk[i].rn[1],&disk[i].rn[2],&disk[i].u[0],&disk[i].u[1],&disk[i].u[2]&disk[i].f[0],&disk[i].f[1],&disk[i].f[2]&disk[i].v[0],&disk[i].v[1],&disk[i].v[2]);
        fgets(tt,2000,dataFile);
    }
}

void populateMatrix(int chain) {
    int i,m;
    double val;
    for (i=0;i<numDisks-1;i++) {
        m = chain*numDisks+i;
        val = calc_cos(m);
        Hamiltonian(i,i+1) = val;
        Hamiltonian(i+1,i) = val;
    }
}

double calc_cos(int i) {
    double value;
    vector n, v1, v2; //Normal vector n defined by the vector between disks, v1 and v2 will be disk 1 and 2's normal (f) vector projected onto plane normal to n from which dihedral angle can be calculated
    for (j=0;j<3;j++) {
        n[j] = disk[i+1].rn[j] - disk[i].rn[j];
        v1[j] = disk[i+1].f[j];
        v2[j] = disk[i].f[j];
    }
    vec_norm(n);
    proj_vector(v1,n);
    proj_vector(v2,n);
    value = vec_dot(v1,v2);
    if (value > 1.0) {
        value = 1.0;
    }
    if (value < -1.0) {
        value = -1.0;
    }
    return value;
}

void Sim::proj_vector(vector &a, vector n) {//Project vector a onto the plane described by unit normal vector n, then normalize the new a
    double dot;
    int i;
    vector n2;
    for (i=0;i<3;i++) {
        n2[i] = n[i];
    }
    dot = vec_dot(a,n2);
    for (i=0;i<3;i++) {
        a[i] = a[i] - dot*n2[i];
    }
    vec_norm(a);
}

void populateHistogram() {

}

int main() {
    int i,j;
    init_data();
    init_matrix();

    for (i=0;i<numSnapshots;i++) {
        readSnapshot();
        for (j=0;j<numChains;j++) {
            populateMatrix(j);
            diagonalize();
            populateHistogram();
        }
    }

    normalizeHistogram();
    kill_everything();
    return 0;
}
