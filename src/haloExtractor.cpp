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
/*
  if ( !( (userInput.getCatType()).compare( "short" ) == 0 ) &&
       !( (userInput.getCatType()).compare( "BMD"   ) == 0 ) ){

    std::cout << " Catalog type unsupported: " << userInput.getCatType() << std::endl;
    exit(1);
  }
*/

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
    printf("\n Number of valid halos: %lu\n",N_halos);
  }

  haloInfo *myHalos = new haloInfo[N_halos];
    printf(" Allocated %lu halos\n",N_halos);

  {
    unsigned long old_N_halos = N_halos;
    N_halos = readCatalog( myHalos, &userInput, N_halos );
    if ( old_N_halos != N_halos ){
      printf("\nError: Read mismatch\nN_halos allocated: %lu\nN_halos read in: %lu\n\n", old_N_halos, N_halos);
      exit(1);
    }
  }
  userInput.setNumHalos( N_halos );


  printf("\n Catalog read in complete\n\n\n");


  /////////////////////////////////////////
  ///////////Read in particles/////////////
  /////////////////////////////////////////

  //setPartFile will locate the particle files
  // if the file exists, or create the file
  // based on the PMSS file it thinks we should use
  //If the read in fails, will return a 0 for num particles

  std::cout << " Attempting to read particle file..." << std::endl;


  // Attempts to locate the particle file we will use,
  //  and reads the number of particles to allocate
  long long numParticles = 0;
  numParticles = setPartFile( userInput );


  // If read 0 particles, error
  if ( numParticles == 0 ){

    printf(" Error reading particle file:%s\n", (userInput.getInputPart()).c_str() );
    printf("Specify the file using the variable \"partFile\", using the variable \"snapNum\", and making sure you are in the proper directory\n\n");
    exit(1);

  }


  // Saves this number of particles
  userInput.setNumParticles( numParticles );
  printf("\n Number of particles: %lli\n", userInput.getNumParticles() );


  // Allocate the partcles
  particlePosition *particle = new particlePosition[numParticles];
  printf(" Allocated %lli particles\n  Reading particle file...\n\n", numParticles );


  // Read the particle positions
  numParticles = readParticle( userInput, particle );
  if ( numParticles != userInput.getNumParticles() ){
    printf("Error: mismatch in particle read in\n Allocated: %lli\n Read in  : %lli\n", userInput.getNumParticles(), numParticles );
    exit(1);
  }


  std::cout << " Particle read in complete\n\n" << std::endl;

  std::cout << " Reading in header file..."     << std::endl;
  readHeader  ( userInput );
  std::cout << " Done."                         << std::endl << std::endl;


  /////////////////////////////////////////
  ////////Generate link lists//////////////
  /////////////////////////////////////////


  std::cout << " Setting up cells..." << std::endl;

  //Determines range of diff coordinates
  //setMinMaxParticles( particle, userInput );

  if ( !userInput.setBox( userInput.getFOV() / 2.0 ) ){
  std::cout << " Error setting boxes "  << std::endl;
    exit(1);
  }
  std::cout << " Done."     << std::endl  << std::endl;


  std::cout << " Generating link list..." << std::endl;


/*
std::cout << "Halo  memory: " <<  sizeof( haloInfo         ) * userInput.getNumHalos()     / (1e6) << " Mb" << std::endl;
std::cout << "Part  memory: " <<  sizeof( particlePosition ) * userInput.getNumParticles() / (1e9) << " Gb" << std::endl;
std::cout << "Link  memory: " <<  sizeof( long long        ) * userInput.getNumParticles() / (1e9) << " Gb" << std::endl;
std::cout << "Label memory: " <<  sizeof( long long        ) * userInput.getNtotCell    () / (1e6) << " Mb" << std::endl;
std::cout << "Total memory: " << (sizeof( long long        ) * userInput.getNtotCell    () +
                                  sizeof( long long        ) * userInput.getNumParticles() +
                                  sizeof( particlePosition ) * userInput.getNumParticles() +
                                  sizeof( haloInfo         ) * userInput.getNumHalos()   ) / (1e9) << " Gb" << std::endl;
//*/

  //Make link list between particles
  long long  *linkList = new long long [ userInput.getNumParticles() ];
  long long *labelList = new long long [ userInput.getNtotCell()     ];

  printf("  Allocated link  list of %lli elements\n", userInput.getNumParticles() );
  printf("  Allocated label list of %i elements\n"  , userInput.getNtotCell()     );


  makeLinkList( userInput, particle, linkList, labelList );
  std::cout << " Done."     << std::endl  << std::endl;


  std::cout << " Refining halo list to those within our cells..." << std::endl;
  //Find the halos who's FOV is entirely inside the box
  {
  //Dummy array, useless
  haloInfo duhalos[2];

  N_halos = findBoxHalos( userInput, myHalos, duhalos, 0 ); //Returns the number of halos in the box on the first run through
  }


  haloInfo halos[N_halos];                                  //  Allocates our new halo array
  N_halos = findBoxHalos( userInput, myHalos,   halos, 1 ); //   Fills in our new halo array
  delete [] myHalos;                                        //        Deletes old halo array
  std::cout << "  Allocated " << N_halos    << " halos" << std::endl;

  std::cout << " Done."       << std::endl  << std::endl;


/*
std::cout << "Halo memory: " << sizeof( haloInfo         ) * userInput.getNumHalos()     / (1e3) << " Kb" << std::endl;
std::cout << "Part memory: " << sizeof( particlePosition ) * userInput.getNumParticles() / (1e6) << " Mb" << std::endl;
//*/


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
