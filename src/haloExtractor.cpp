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

  long long numParticles = 0;

  numParticles = setPartFile( userInput );

  printf("\n");
  if ( numParticles == 0 ){

    printf(" Error reading particle file:%s\n", (userInput.getInputPart()).c_str() );
    printf("Specify the file using the variable \"partFile\", using the variable \"snapNum\", and making sure you are in the proper directory\n\n");
    exit(1);

  }
  userInput.setNumParticles( numParticles );
  //Once particles are confirmed, read into this C side
  //particlePosition particle[numParticles];                       //This method segfaults due to compiler limit

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


//Loop over halos, find boxes to check with
//Z axis should check against long integration lengths
//Loop once to find, generate link lists for each set
//sphere, box, integration lengths (if applicable)
//Integration length images first, write file
//Box, write file
//Rotate sphere for different views, write file


//For now do arrays, figure out the fits writing later

//32 boxes
  linkHaloParticles( userInput, halos, particle, labelList, linkList );
/*
int Nlx = userInput.getNlx(); //Minimum box numbers (should be 0)
int Nly = userInput.getNly();
int Nlz = userInput.getNlz();
int Nrx = userInput.getNrx(); //Maximum box numbers
int Nry = userInput.getNry();
int Nrz = userInput.getNrz();

float xMin = userInput.getXmin();
float yMin = userInput.getYmin();
float zMin = userInput.getZmin();
float cell = userInput.getCell();
float FOV  = userInput.getFOV() ;

for ( int i = 0; i < userInput.getNumHalos(); ++i ){

  float x = halos[i].getX();
  float y = halos[i].getY();
  float z = halos[i].getZ();

  //Box to the left and right of the box the center of the halo is located in
  int xLeft  = std::min( std::max( Nlx, int( ( halos[i].getX() - xMin ) / cell - 1 ) ), Nrx );
  int yLeft  = std::min( std::max( Nly, int( ( halos[i].getY() - yMin ) / cell - 1 ) ), Nry );
  int zLeft  = std::min( std::max( Nlz, int( ( halos[i].getZ() - zMin ) / cell - 1 ) ), Nrz );
  int xRight = std::min( std::max( Nlx, int( ( halos[i].getX() - xMin ) / cell + 1 ) ), Nrx );
  int yRight = std::min( std::max( Nly, int( ( halos[i].getY() - yMin ) / cell + 1 ) ), Nry );
  int zRight = std::min( std::max( Nlz, int( ( halos[i].getZ() - zMin ) / cell + 1 ) ), Nrz );


  //Number of particles in each set
  unsigned long N_sphere(0);
  unsigned long N_box   (0);
  unsigned long N_i1    (0);
  unsigned long N_i2    (0);
  unsigned long N_i3    (0);


  //Loop over cells possibly containing particles of our halo, and count
  // the number that belong in each particle set
  for ( int xIndex = xLeft; xIndex <= xRight; ++xIndex ){
  for ( int yIndex = yLeft; yIndex <= yRight; ++yIndex ){
  for ( int zIndex = zLeft; zIndex <= zRight; ++zIndex ){

    long cellIndex = xIndex +                          //Index for the cell for labellist
                     yIndex * Nrx +
                     zIndex * Nrx * Nry;

    long long particleIndex = labelList[ cellIndex ];  //The index for the first particle in the link list

    //Loop over link list
    while ( particleIndex != -1 ){

      float partX = particle[ particleIndex ].x_pos;
      float partY = particle[ particleIndex ].y_pos;
      float partZ = particle[ particleIndex ].z_pos;

      //Check if in sphere, Rm > sqrt( ( x_h-x_p )^2 + ...
      if ( ( halos[i].getRm() * halos[i].getRm() ) > ( ( halos[i].getX() - partX ) * ( halos[i].getX() - partX ) +
                                                       ( halos[i].getY() - partY ) * ( halos[i].getY() - partY ) +
                                                       ( halos[i].getZ() - partZ ) * ( halos[i].getZ() - partZ ) ) ){
        ++N_sphere;
      }

      //Particle in the box frame
      else if ( ( (halos[i].getX()-FOV) < partX ) && ( partX < (halos[i].getX()+FOV) ) &&
                ( (halos[i].getY()-FOV) < partY ) && ( partY < (halos[i].getY()+FOV) ) &&
                ( (halos[i].getZ()-FOV) < partZ ) && ( partZ < (halos[i].getZ()+FOV) ) ){
        ++N_box;
      }
//Check if in box
//Integration length, needs modification of z index


      particleIndex = linkList[ particleIndex ];
    }
  }
  }
  }
N_box = 1;


  //Only continue if we found particles for our halo
  if ( ( N_sphere > 0 ) && ( N_box > 0 ) ){


N_box = 0;
    //Need to allocate for indexes
    long long *sphereIndexes = new long long [ N_sphere ];
    long long *   boxIndexes = new long long [ N_box    ];
    long long *    i1Indexes;
    long long *    i2Indexes;
    long long *    i3Indexes;

    //If we were able to complete any integration, use it
    if ( N_i1 > 0 ) i1Indexes = new long long [ N_i1 ];
    if ( N_i2 > 0 ) i2Indexes = new long long [ N_i2 ];
    if ( N_i3 > 0 ) i3Indexes = new long long [ N_i3 ];

    long long  sCounter(0);
    long long  bCounter(0);
    long long i1Counter(0);
    long long i2Counter(0);
    long long i3Counter(0);

    //Loop over cells again to find indexes of particles
    for ( int xIndex = xLeft; xIndex <= xRight; ++xIndex ){
    for ( int yIndex = yLeft; yIndex <= yRight; ++yIndex ){
    for ( int zIndex = zLeft; zIndex <= zRight; ++zIndex ){

      long cellIndex = xIndex +                          //Index for the cell for labellist
                       yIndex * Nrx +
                       zIndex * Nrx * Nry;

      long long particleIndex = labelList[ cellIndex ];  //The index for the first particle in the link list

      //Loop over link list
      while ( particleIndex != -1 ){

        float partX = particle[ particleIndex ].x_pos;
        float partY = particle[ particleIndex ].y_pos;
        float partZ = particle[ particleIndex ].z_pos;

        //Check if in sphere, Rm > sqrt( ( x_h-x_p )^2 + ...
        if ( ( halos[i].getRm() * halos[i].getRm() ) > ( ( halos[i].getX() - partX ) * ( halos[i].getX() - partX ) +
                                                         ( halos[i].getY() - partY ) * ( halos[i].getY() - partY ) +
                                                         ( halos[i].getZ() - partZ ) * ( halos[i].getZ() - partZ ) ) ){
          sphereIndexes[ sCounter ] = particleIndex;
          ++sCounter;
        }

        //Particle in the box frame
        else if ( ( (halos[i].getX()-FOV) < partX ) && ( partX < (halos[i].getX()+FOV) ) &&
                  ( (halos[i].getY()-FOV) < partY ) && ( partY < (halos[i].getY()+FOV) ) &&
                  ( (halos[i].getZ()-FOV) < partZ ) && ( partZ < (halos[i].getZ()+FOV) ) ){
          boxIndexes[ bCounter ] = particleIndex;
          ++bCounter;
        }


      particleIndex = linkList[ particleIndex ];
    }

    }
    }
    }




  }
}

//*/







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
