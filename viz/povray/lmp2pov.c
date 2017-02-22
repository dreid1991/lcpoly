// convert coordinates and quaternions from LAMMPS format to coordinates and vectors in format understood by POV-ray
//   read <'oneframe'.lammpstrj>
//   convert quaternions to vectors
//   write coordinates and anlges to use in POV-ray (additional macros for ellipsoid/spherocylinder needs to be defined in <head.pov>)

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#define  PI 3.14159265358979323846

void read_lammpstrj(char *lmp);
void quat2vec();
void write_pov(char *pov);
//int nrsr_int(double argm);
//void EndLine(FILE *file);

struct atom
{
    int id, type;
    double x, y, z;
    double qw, qi, qj, qk;
} *cord;

struct rot
{
    double x, y, z;
} *vec;

int natom, nsub;
double boxs, xlo, xhi, ylo, yhi, zlo, zhi, Sx, Sy, Sz;
char style;

int main(int argc, char *argv[])
{
   // check usage:
    if(argc!=6)
    {
        printf("usage: %s <frame.lammpstrj>  Sx Sy Sz E|C|P \n", argv[0]);
        printf("       reads one frame of GB ellipsoids in *.lammpstrj format\n");
        printf("       writes converted <out.pov> for use with POV-ray\n");
        printf("       Sx, Sy, Sz - tree axes of ellopsoid defining its shape\n");
        printf("       E|C|P - representation style for ellipsoid|spherocylinder|spheroplatelet\n");
        exit(0);
    }

   // read one frame:
    read_lammpstrj(argv[1]);

   // set the shape of ellipsoids:
    Sx = atof(argv[2]);
    Sy = atof(argv[3]);
    Sz = atof(argv[4]);
 
   // representation style:
    style = argv[5][0];
    switch(style) {
     // as ellipsoids: 
      case 'E':
        printf("E - representation style is set to Ellipsoid\n");
        break;

     // as spherocylinders: (assuming z- or x-axis as longest)
      case 'C':
        printf("C - representation style is set to sphero-Cylinder\n");
        break;

     // as spheroplatelet: (assuming x- and y-axes are long)
      case 'P':
        printf("P - representation style is set to sphero-Platelet\n");
        break;
      default :
        printf("Invalid option for represetation style. It can be only E|C|P\n");
        exit(0);
    }

   // convert quaternions to vectors:
    //quat2vec();

   // write ouput:
    write_pov("out.pov");
}

void read_lammpstrj(char *lmp)
{
    int i, id, type;
    double x, y, z, ax,ay,az;
    char c, text[100];
    FILE *input;

    input = fopen(lmp, "r"); // open input to read
    printf("reading input file <%s>\n", lmp);

    for(i=0; i<3; i++) fgets(text, 100, input); // skip lines

    fscanf(input, "%d", &natom);
    printf("number of atoms: %d\n", natom);

    for(i=0; i<2; i++) fgets(text, 100, input); // skip lines again

   // simulations box sizes:
    fscanf(input, "%lf %lf", &xlo, &xhi);
    fscanf(input, "%lf %lf", &ylo, &yhi);
    fscanf(input, "%lf %lf", &zlo, &zhi);
    printf("box dimentions are:\n  xlo xhi   %lf  %lf\n  ylo yhi   %lf  %lf\n  zlo zhi   %lf  %lf\n", xlo, xhi, ylo, yhi, zlo, zhi);

    for(i=0; i<2; i++) fgets(text, 100, input); // skip lines

   // allocate memory for coordinates:
    cord = calloc(natom, sizeof(struct atom));
    vec = calloc(natom, sizeof(struct rot));

   // read coordinates and quaternions:
    for(i=0; i<natom; i++)
    {
        fscanf(input, "%d %d %lf %lf %lf %lf %lf %lf", &id, &type, &x, &y, &z, &ax, &ay, &az);
        cord[i].id = id;
        cord[i].type = type;
        cord[i].x = x;
        cord[i].y = y;
        cord[i].z = z;
        vec[i].x = ax;
        vec[i].y = ay;
        vec[i].z = az;
        /*cord[i].qw = qw;
        cord[i].qi = qi;
        cord[i].qj = qj;
        cord[i].qk = qk;*/
    }
    fgets(text, 100, input);

    fclose(input);
}

// convert quaternions to vectors:
void quat2vec()
// the 3x3 rotation matrix is expressed in terms of quaternion components {qw,qi,qj,qk} as:
//   
//   |1-2*(qj*qj+qk*qk)    2*(qi*qj-qk*qw)    2*(qi*qk+qj*qw) |
// R=|  2*(qi*qj+qk*qw)  1-2*(qi*qi+qk*qk)    2*(qj*qk-qi*qw) |
//   |  2*(qi*qk-qj*qw)    2*(qj*qk+qi*qw)  1-2*(qi*qi+qj*qj) |
//
// the column-vertor a={ax;ay;az} needs to be multiplied by matrix R from left side: Ra = b
{
    int i;
    double qw, qi, qj, qk;

   // allocate memory for vectors:
    vec = calloc(natom, sizeof(struct rot));

   // rotating vector {1,0,0}  (must exclude degeneracy in describing ellipsoid. {1,0,0} is not good for 1-1-3) 
   //    in principle is shouldn't matter what definition in LAMMPS was used: 3-1-1 or 1-1-3
   //    the sphere will be stretched to ellipsoid accordingly (along x- or z-axis) (THAT IS TO BE FIGURED OUT)
    for(i=0;i<natom;i++)
    {
        vec[i].x = 1 - 2*(pow(cord[i].qj,2) + pow(cord[i].qk,2));
        vec[i].y = 2*(cord[i].qi*cord[i].qj + cord[i].qk*cord[i].qw);
        vec[i].z = 2*(cord[i].qi*cord[i].qk - cord[i].qj*cord[i].qw);
//        qw = cord[i].qw; 
//        qi = cord[i].qi;
//        qj = cord[i].qj;
//        qk = cord[i].qk;
//        vec[i].x = 1 - 2*(qj*qj + qk*qk)  +  2*(qi*qj - qk*qw)  +  2*(qi*qk + qj*qw);
//        vec[i].y = 2*(qi*qj + qk*qw)  +  1 - 2*(qi*qi + qk*qk)  +  2*(qj*qk - qi*qw);
//        vec[i].z = 2*(qi*qk - qj*qw)  +  2*(qj*qk + qi*qw)  +  1 - 2*(qi*qi + qj*qj);
    }

    printf("quaternions have been converted to vectors\n");
}

// write output in POV-Ray format:
void write_pov(char *pov)
// output format is for macros that is to be defined later in other script (head.pov):
//    ellipsoid(  Position,        Rotation,           Color,       Scale)
//               <x, y, z> , <vec_x, vec_y, vec_z>,  <r, g, b>,  <Sx, Sy, Sz>
{
    int i;
    double x, y, z, ax, ay, az, r, g, b, norm;
    FILE *output;

    output = fopen(pov, "w"); // open output file 

    for(i=0; i<natom; i++)
    {
        x = cord[i].x;
        y = cord[i].y;
        z = cord[i].z;
       // vectors:
        ax = vec[i].x;
        ay = vec[i].y;
        az = vec[i].z;
        norm = sqrt(ax*ax + ay*ay + az*az);
       // color scheme:
        r = fabs(ax/norm);
        g = fabs(ay/norm);
        b = fabs(az/norm);
       // print to file:
        if(cord[i].type == 1)
        {
            switch(style) {
          // as ellipsoids: 
              case 'E':
                fprintf(output, "ellipsoid(<%lf,%lf,%lf>,<%lf,%lf,%lf>,<%lf,%lf,%lf>,<%lf,%lf,%lf>)\n", x, y, z, ax, ay, az, r, g, b, Sx, Sy, Sz);
                break;

          // as spherocylinders: (assuming z- or x-axis as longest)
              case 'C':
                fprintf(output, "spherocylinder(<%lf,%lf,%lf>,<%lf,%lf,%lf>,<%lf,%lf,%lf>,0.5,%lf)\n", x, y, z, ax, ay, az, r, g, b, Sx-1);
                break;

          // as spheroplatelet: (assuming x- and y-axes are long)
              case 'P':
                fprintf(output, "spheroplatelet(<%lf,%lf,%lf>,<%lf,%lf,%lf>,<%lf,%lf,%lf>,%lf,%lf)\n", x, y, z, ax, ay, az, r, g, b, 0.5*(Sx-Sz), 0.5*Sz);
                break;
            }
        }

        // properties for substrate atoms:
        else fprintf(output, "ellipsoid(<%lf,%lf,%lf>,<%lf,%lf,%lf>,<%lf,%lf,%lf>,<%lf,%lf,%lf>)\n", x, y, z, ax, ay, az, r, g, b, 1.0, 1.0, 1.0);
    }
    printf("output is written to <%s>\n", pov);

    fclose(output); 
}
