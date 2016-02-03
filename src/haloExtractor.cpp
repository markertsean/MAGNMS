#include <iostream>
#include <fstream>
#include <cstring>
#include <sys/stat.h>


#include "halo_extraction_classes.h"
#include "read_files.h"
#include "link_halos.h"

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

//  haloInfo myHalos[N_halos];
  haloInfo *myHalos = new haloInfo[N_halos];

  {
    unsigned long old_N_halos = N_halos;
    N_halos = readCatalog( myHalos, userInput, N_halos );
    if ( old_N_halos != N_halos ){
      printf("\nError: Read mismatch\nN_halos allocated: %lu\nN_halos read in: %lu\n\n", old_N_halos, N_halos);
      exit(1);
    }
  }
  userInput.setNumHalos( N_halos );


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

  //setPartFile will locate the particle files
  // if the file exists, or create the file
  // based on the PMSS file it thinks we should use
  //If the read in fails, will return a 0 for num particles

  unsigned long numParticles = 0;

  numParticles = setPartFile( userInput );

  printf("\n");

  if ( numParticles == 0 ){

    printf(" Error reading particle file:%s\n", (userInput.getInputPart()).c_str() );
    printf("Specify the file using the variable \"partFile\", using the variable \"snapNum\", and making sure you are in the proper directory\n\n");
    exit(1);

  }
  userInput.setNumParticles( numParticles );

  //Once particles are confirmed, read into this C side
  particlePosition particle[numParticles];


  readParticle( userInput, particle );

  setMinMaxParticles( particle, userInput );


  /////////////////////////////////////////
  ////////Write Halo Particle Files////////
  /////////////////////////////////////////

//Shitton of particles
//Loop over halos, need some sort of flag to mark as used/unused
//Check particles for boundaries in x, y, z, will vary on read in
//  reduce halos list to those within boundary
//float particlePosition::x_min=0;

  haloInfo *halos;
  {
    //Array to track how many halos are valid within the box
    short inBox[ userInput.getNumHalos() ];
    unsigned long sum(0);
std::cout <<"166"<<std::endl;

    //Cycle through halos to determine if in our box
    for ( int i  = 0 ; i < userInput.getNumHalos(); ++i ){
        inBox[i] = 0;

      //If the halo is within the coordinates we have, and the FOV doesn't leave the box, save it
      if ( ( userInput.getXmin() < myHalos[i].getX() - userInput.getFOV() ) && ( myHalos[i].getX() + userInput.getFOV() < userInput.getXmax() ) &&
           ( userInput.getYmin() < myHalos[i].getY() - userInput.getFOV() ) && ( myHalos[i].getY() + userInput.getFOV() < userInput.getYmax() ) &&
           ( userInput.getZmin() < myHalos[i].getZ() - userInput.getFOV() ) && ( myHalos[i].getZ() + userInput.getFOV() < userInput.getZmax() )   ){

        inBox[i] = 1;
      }

      //Sum is the number of halos in the box
      sum += inBox[i];
    }
std::cout <<"183"<<std::endl;

    //Allocate the new halo list
    halos = new haloInfo[  sum ];

    //Go through and copy over the halos in the box to the new halo list
    unsigned long counter(0);
    for ( int i = 0 ; i < userInput.getNumHalos(); ++i ){
      if ( inBox[i] == 1 ){
        halos[ counter ] = myHalos[ i ];
        ++counter;
std::cout << myHalos[i].getX() << std::endl;
      }
    }
std::cout <<"197"<<sum<<std::endl;

    userInput.setNumHalos( sum );
  }

  delete [] myHalos;

for ( int i = 0; i < userInput.getNumHalos(); ++i ){
std::cout << halos[i].getX() << std::endl;
}



//Make link list
int blah = 8;
  unsigned long LinkList[ userInput.getNumParticles() ];
//  int label[cellnum,cellnum,cellnum]




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



  delete [] halos;
}
