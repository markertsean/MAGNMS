#include <iostream>
#include <cstring>
#include <sys/stat.h>

#include "halo_extraction_classes.h"
#include "read_files.h"

//Read halo catalog
//Select halos within criterion
//Write useable catalog, if one doesn't exist
//Create directories for infofile, particle positions
//Go through particle catalog, extracting

//Need read functions for particle file, halo catalog
//Identify catalog/particle file type

//Function to to create directories using file names




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
    printf("  Halo Catalog       : %s\n  Particle File      : %s\n  Header Directory   : %s\n  Particle Directory : %s\n\n",
    (userInput.getInputCatalog()).c_str(),
    (userInput.getInputPart   ()).c_str(),
    (userInput.getHeaderDir   ()).c_str(),
    (userInput.getParticleDir ()).c_str());
  }
  else {
    printf("Error in reading input file, aborting\n\n");
    exit(1);
  }


  ///////////////////////////////////////////
  ////////Read in the catalog file///////////
  ///////////////////////////////////////////


  //Read the number of valid halos, for allocation
  unsigned long N_halos =0;
  {
    haloInfo myHalos[2];
    N_halos = readCatalog( myHalos, userInput, N_halos );
    printf(" Number of halos read: %lu\n",N_halos);
  }

  haloInfo myHalos[N_halos];
  {
    unsigned long old_N_halos = N_halos;
    N_halos = readCatalog( myHalos, userInput, N_halos );
    if ( old_N_halos != N_halos ){
      printf("\nError: Read mismatch\nN_halos allocated: %lu\nN_halos read in: %lu\n\n", old_N_halos, N_halos);
      exit(1);
    }
  }

  printf("\n Catalog read in complete\n\n");


  /////////////////////////////////////////
  //////Write header/part directory////////
  /////////////////////////////////////////

  {
    struct stat sb;
    int maxStrLength=400;
    char str[maxStrLength];

    //If header directory exists, we use it
    //If not, create it
    sprintf( str, "mkdir %s", (userInput.getHeaderDir()).c_str() );
    if ( stat((userInput.getHeaderDir()).c_str(), &sb) == 0 ){
      std::cout << " Found directory: " << userInput.getHeaderDir() << std::endl;;
    } else {
      std::cout << " Writing directory: " << userInput.getHeaderDir() << std::endl;
      system(str);
    }

    sprintf( str, "mkdir %s", (userInput.getParticleDir()).c_str() );
    if ( stat((userInput.getParticleDir()).c_str(), &sb) == 0 ){
      std::cout << " Found directory: " << userInput.getParticleDir() << std::endl;;
    } else {
      std::cout << " Writing directory: " << userInput.getParticleDir() << std::endl;
      system(str);
    }

    printf("\n\n");
  }


  /////////////////////////////////////////
  ///////////Read in particles/////////////
  /////////////////////////////////////////







  /////////////////////////////////////////
  //////////Write header files/////////////
  /////////////////////////////////////////

  //Catalog
  //ID number
  //Halo position
  //Particle number from catalog
  //Mass of particles
  //Halo mass
  //Rmax
  //Concentration
  //b/a ratio
  //c/a ratio
  //
  //Rotation vector for along LOS
  //Rotation vector for perp to LOS
  //Number particles using
  //Number particles in box
  //Whether could create LOS cones
  //Largest LOS cones

  //Have input files and directories to write to.
  //Need to:
  //        Write overarching header file, with used halos
  //        Read in particles, maybe use link list to speed up process?




}
