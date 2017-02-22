################################################################################
# this makefile will prepare and render POV-Ray scene of Gay-Berne ellipsoids 
################################################################################

### TO DO later: 
#  1) interprete optional input parameters not as INPUT=file.name,
#    but as 1st, 2nd, 3rd arguments
#  2) print usage message (???)
#        make   - will use default values for parameters
#        make  frame.lammpstrj  3  1  1  low   - list of parameters
#                 INPUT        SX SY SZ  RES
################################################################

# filename for one frame in LAMMPS format:
INPUT = disks.xyz 

# default shape of ellipoidal particles:
SX = 0.5
SY = 1
SZ = 1
# for other shapes either correct values above or use make with options:
#   [/workdir]$  make  SX=1 SY=1 SZ=3

# representation style E|C|P (ellipsoid, spherocylinder, spheroplatelet):
STYLE=E

# default resolution for output image:
RES = low

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
all: render

render: scene.pov RES.ini
	povray RES[$(RES)] scene.pov

# prepare scene configuration for POV-ray: 
scene.pov: head.pov out.pov
	cat head.pov out.pov > scene.pov

# convert frame of GB particle in LAMMPS format to POV-ray (quaternions to angles etc. :
out.pov: lmp2pov
	./lmp2pov $(INPUT) $(SX) $(SY) $(SZ) $(STYLE)

# compile c-code - routing for format conversion:
lmp2pov: lmp2pov.c
	gcc -o lmp2pov lmp2pov.c -lm

#----------------------------------------
clean:
	rm out.pov scene.pov

clean-all:
	rm out.pov scene.pov scene.png lmp2pov 
