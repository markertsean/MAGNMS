/***************************************************************************************

This program takes an input file (extractInfo.dat or command line argument) for input
parameters to write FITS mass maps from simulations. The process is:

1. Read the user input
2. Perform internal checks for directories
3. Read the halo catalog, write a short file of valid halos for faster read in
4. Read in the particles, save as array into memory
5. Generate link list of particles
6. Link particles and halos
   -Particles go in Rvir, FOV box, or integration list sets
   -Writes images with headers and FITS mass maps, no smoothing
   -Currently does nothing with triaxiality

The last step tries to avoid allocating as much memory as possible, due to the extreme number of particles.


****************************************************************************************/

#include <iostream>
#include <fstream>
#include <cstring>
#include <sys/stat.h>

#include <CCfits/CCfits>
#include <cmath>

#include "halo_extraction_classes.h"
#include "read_files.h"
#include "link_halos.h"


unsigned long inputInfo::Num_files = 0;


int main( int arg, char ** argv ){

  //////////////////////////////////////
  //////////User info read in///////////
  //////////////////////////////////////


  //Stores the user input from the input file
  inputInfo userInput;

  //User can specify input file via command line, default: extractInfo.dat
  if ( arg > 1 ){
    userInput.setReadFile( argv[1] );
  }
  printf("\n Using input file: %s\n", (userInput.getReadFile()).c_str());

  //Attempts to read the input file, and displays the resultant info
  if ( readUserInput( userInput.getReadFile(), userInput ) ){
    printf("  Halo Catalog       : %s\n  Particle File      : %s\n  Mass map Directory : %s\n\n",
    (userInput.getInputCatalog()).c_str(),
    (userInput.getInputPart   ()).c_str(),
//    (userInput.getHeaderDir   ()).c_str(),
    (userInput.getParticleDir ()).c_str());
  }
  else {
    printf("Error in reading input file, aborting\n\n");
    exit(1);
  }

  if ( !( (userInput.getCatType()).compare( "short" ) == 0 ) &&
       !( (userInput.getCatType()).compare( "BMD"   ) == 0 ) ){

    std::cout << " Catalog type unsupported: " << userInput.getCatType() << std::endl;
    exit(1);
  }


  /////////////////////////////////////////
  //////Write header/part directory////////
  /////////////////////////////////////////

  {
    struct stat sb;
    int maxStrLength=400;
    char str[maxStrLength];

    /*
    //If header directory exists, we use it
    //If not, create it
    sprintf( str, "mkdir %s", (userInput.getHeaderDir()).c_str() );
    if ( stat((userInput.getHeaderDir()).c_str(), &sb) == 0 ){
      std::cout << " Found directory: " << userInput.getHeaderDir() << std::endl;;
    } else {
      std::cout << " Writing directory: " << userInput.getHeaderDir() << std::endl;
      system(str);
    }
    */

    sprintf( str, "mkdir %s", (userInput.getParticleDir()).c_str() );
    if ( stat((userInput.getParticleDir()).c_str(), &sb) == 0 ){
      std::cout << "  Found directory    : " << userInput.getParticleDir() << std::endl;;
    } else {
      std::cout << "  Writing directory  : " << userInput.getParticleDir() << std::endl;
      system(str);
    }

    printf("\n\n");
  }



  ///////////////////////////////////////////
  ////////Read in the catalog file///////////
  ///////////////////////////////////////////

      std::cout << " Attempting to read halo catalog..." << std::endl;


  //Read the number of valid halos, for allocation
  unsigned long N_halos =0;
  {
    haloInfo myHalos[2];
    N_halos = readCatalog( myHalos, &userInput, N_halos );
    printf("\n   Number of valid halos: %lu\n\n",N_halos);
  }

  haloInfo *myHalos = new haloInfo[N_halos];

  {
    unsigned long old_N_halos = N_halos;
    N_halos = readCatalog( myHalos, &userInput, N_halos );
    if ( old_N_halos != N_halos ){
      printf("\nError: Read mismatch\nN_halos allocated: %lu\nN_halos read in: %lu\n\n", old_N_halos, N_halos);
      exit(1);
    }
  }
  userInput.setNumHalos( N_halos );


  printf("\n Catalog read in complete\n\n");


  /////////////////////////////////////////
  ///////////Read in particles/////////////
  /////////////////////////////////////////

  //setPartFile will locate the particle files
  // if the file exists, or create the file
  // based on the PMSS file it thinks we should use
  //If the read in fails, will return a 0 for num particles

  long long numParticles = 0;

  numParticles = setPartFile( userInput );

  printf("\n");
  if ( numParticles == 0 ){

    printf(" Error reading particle file:%s\n", (userInput.getInputPart()).c_str() );
    printf("Specify the file using the variable \"partFile\", using the variable \"snapNum\", and making sure you are in the proper directory\n\n");
    exit(1);

  }
  userInput.setNumParticles( numParticles );


  particlePosition *particle = new particlePosition[numParticles]; //This does not segfault


  readParticle( userInput, particle );


  /////////////////////////////////////////
  ////////Generate link lists//////////////
  /////////////////////////////////////////

  std::cout << " Setting up cells..." << std::endl;

  //Determines range of diff coordinates
  setMinMaxParticles( particle, userInput );

  //Find the halos who's FOV is entirely inside the box
  {
  //Dummy array, useless
  haloInfo duhalos[2];

  N_halos = findBoxHalos( userInput, myHalos, duhalos, 0 ); //Returns the number of halos in the box
  }

  haloInfo halos[N_halos];                                  //Allocates our new halo array
  N_halos = findBoxHalos( userInput, myHalos,   halos, 1 ); //Fills in our new halo array
  delete [] myHalos;                                        //Deletes old halo array


  if ( !userInput.setBox( userInput.getFOV() / 2.0 ) ){
    std::cout << " Error setting boxes "  << std::endl;
    exit(1);
  }

  std::cout << " Done."     << std::endl  << std::endl;


  std::cout << " Generating link list..." << std::endl;

  //Make link list between particles
  long long  *linkList = new long long [ userInput.getNumParticles() ];
  long long *labelList = new long long [ userInput.getNtotCell()     ];

  makeLinkList( userInput, particle, linkList, labelList );
  std::cout << " Done."     << std::endl  << std::endl;


  /////////////////////////////////////////
  ////////Link halos, write files//////////
  /////////////////////////////////////////

  // Goes through link list, checking particles against halos
  // Locates particles in halo's Rvir
  //  box FOV x FOV x FOV
  //  box FOV x FOV x integration lengths
  // Writes image to FITS file

  std::cout << " Generating FITS images..." << std::endl << std::endl;
  linkHaloParticles( userInput, halos, particle, labelList, linkList );
  std::cout << " Done."                     << std::endl << std::endl;

  std::cout << "Wrote " << userInput.getNumFiles() << " Files " << std::endl;

}
