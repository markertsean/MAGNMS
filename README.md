The Mass mAp Generator for the Multidark Suite is designed to take a PMss file from
the MultiDark suite, and convert it to a FITS surface density map compatible with 
the GLAMER code (http://glenco.github.io/glamer/). 

The code uses CMake to compile the code. A CMakeLists.txt is included, but will 
require modification for each user.

This program takes an input file (extractInfo.dat or file given as command line argument) 
for input parameters to write FITS mass maps from simulations. Each argument must be listed,
and the value seperated by whitespace. Accepted variables:

haloFile        -               Halo catalog to use, must have same format as MDPlanck CatshortW files
partFile        -               PMss particle file to use
headFile        -               Header file to use for PMss file, if seperate from PMss file
partDir         -               Directory of particle file, useful for when multiple PMss particles being read in
headDir         -               Directory of header file, useful for when multiple PMss particles being read in
catType         - Default MD  , CURRENTLY ONLY SUPPORT FOR MDP (MultiDark Planck)
integAxis       - Default z   , support for x, y, and z axis
minMass         - Defualt 1e15, lower halo mass limit to read from file
maxMass         -               Upper halo mass to read from file
rMult           - Default 1.0 , sphere images extend to rMult * r_vir
rConv           - Default 1e-3, conversion factor between r_vir units and coordinates
FOV             - Default 8.0 , field of view in h^-1 Mpc
integLength     - Default 400 , maximum integration length in h^-1 Mpc
integStep       - Default 0.2 , integration length step in dex
useShortCatalog - Default 1   , 0 or 1 to turn on or off reading/writing a short list of particles to speed up read in
snapNum         -               PMss snapshot number to search for
N_pixels_h      - Default 1024, number of horizontal pixels in FITS image
N_pixels_v      - Default 1024, number of vertical pixels in FITS image
physicalSize    -               Misleading title, angular field of view in degrees of image
redshift        -               Redshift of halo, does not need to match catalog
PMssFirst       -               First node of PMss files to search for
PMssLast        -               Last  node of PMss files to search for

The procedure of the code is outlined below:

1. Read the user input from extractInfo.dat

2. Perform internal checks for directories

3. Read the halo catalogs, write a short file of valid halos for faster read in based on user constraints

4. Read in the particles from PMss files, save as array into memory. The PMss files will be searched for 
in the partDir directory, and search based on recognized file patterns with the snapNum and PMss first and last values.
Additionally, the memory usage here is extreme, be mindful of memory usage of your computer before running.

5. Generate link list of particles

6. Link particles and halos, writing output FITS files
   -Particles go in Rvir, FOV box, or integration list sets
   -Writes images with headers and FITS mass maps, with no smoothing
   -Calculates the orientation of the halo, writes to the header
