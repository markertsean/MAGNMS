/***************************************************************************************

This program takes an input file (extractInfo.dat or given as command line argument) for input
parameters to write FITS mass maps from simulations. The process is:

1. Read the user input
2. Perform internal checks for directories
3. Read the halo catalog, write a short file of valid halos for faster read in based on user constraints
4. Read in the particles, save as array into memory
5. Generate link list of particles
6. Link particles and halos, writing output FITS files
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
#include <ctime>

#include "halo_extraction_classes.h"
#include "read_files.h"
#include "link_halos.h"
#include "triaxiality.h"


// Keeps track of the number of FITS files written
unsigned long inputInfo::Num_files = 0;


std::string logFileName = "";


int main( int arg, char ** argv ){


  // Get the time the code started at

  int execution_start = clock();
  time_t      nowTime =      time(        0 );
  tm       *startTime = localtime( &nowTime );

  // Generates the log file name
  {
    char logFileNameC[100];

    struct stat sb;
    char str[] = "mkdir logfiles";

    // Create logfiles directory, if one does not exist
    if ( stat( "logfiles/", &sb) != 0 ){
      system( str );
    }

    // Log file name is mostly the date
    sprintf( logFileNameC, "logfiles/logfile_MAGNMS.%4i.%02i.%02i.%02i.%02i.%02i.log",
      (*startTime).tm_year+1900,
      (*startTime).tm_mon ,
      (*startTime).tm_mday,
      (*startTime).tm_hour,
      (*startTime).tm_min ,
      (*startTime).tm_sec );

    logFileName= std::string(logFileNameC) ;

    // First line of the log file is the time
    logMessage( (std::string(                "Code initialized at ")+
                 std::to_string( (long long) (*startTime).tm_hour  )+
                 std::string(                ":"                   )+
                 std::to_string( (long long) (*startTime).tm_min   )+
                 std::string(                ":"                   )+
                 std::to_string( (long long) (*startTime).tm_sec   )));
  }


  //////////////////////////////////////
  //////////User info read in///////////
  //////////////////////////////////////


  // Stores the user input from the input file
  inputInfo userInput;

  // User can specify input file via command line, default: extractInfo.dat
  if ( arg > 1 ){
    userInput.setReadFile( argv[1] );
  }


  printf("\n Using input file: %s\n", (userInput.getReadFile()).c_str());

  logMessage( std::string("Input file: ") +
                  userInput.getReadFile() );


  // Attempts to read the input file, and displays the resultant info
  if ( readUserInput( userInput.getReadFile(), userInput ) ){

    printf("  Halo Catalog       : %s\n  Particle File      : %s\n  Mass map Directory : %s\n\n",
    (userInput.getInputCatalog()).c_str(),
    (userInput.getInputPart   ()).c_str(),
    (userInput.getParticleDir ()).c_str());

  }
  else {

    printf("Error in reading input file, aborting\n\n");

    logMessage( "Aborting" );

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
//*/

    sprintf( str, "mkdir %s", (userInput.getParticleDir()).c_str() );
    if ( stat((userInput.getParticleDir()).c_str(), &sb) == 0 ){
      std::cout << "  Found directory    : " << userInput.getParticleDir() << std::endl;
    } else {
      std::cout << "  Writing directory  : " << userInput.getParticleDir() << std::endl;
      system(str);
    }

    logMessage( std::string("Using particle directory: ") +
                std::string( userInput.getParticleDir() ));

    printf("\n\n");
  }



  ///////////////////////////////////////////
  ////////Read in the catalog file///////////
  ///////////////////////////////////////////

      std::cout << " Attempting to read halo catalog..." << std::endl;


  // Read the number of valid halos, for allocation
  unsigned long N_halos =0;
  {
    haloInfo myHalos[2];
    N_halos = readCatalog( myHalos, &userInput, N_halos );
    printf("\n Number of valid halos: %lu\n",N_halos);
  }

  // myHalos is an array containing all relevant information about the halos

  haloInfo *myHalos = new haloInfo[N_halos];
    printf(" Allocated %lu halos\n",N_halos);

  logMessage( std::string(            "Allocated ") +
              std::to_string( (long long) N_halos ) +
              std::string(        " halos in main") );


  {

    unsigned long old_N_halos = N_halos;
    N_halos = readCatalog( myHalos, &userInput, N_halos );

    if ( old_N_halos != N_halos ){
      printf("\nError: Read mismatch\nN_halos allocated: %lu\nN_halos read in: %lu\n\n", old_N_halos, N_halos);

      logMessage( std::to_string( (long long) old_N_halos     ) +
                  std::string(            " halos allocated, ") +
                  std::to_string( (long long) N_halos         ) +
                  std::string(  " halos read the second time. Aborting") );

      exit(1);
    }
  }

  userInput.setNumHalos( N_halos );


  printf("\n Catalog read in complete\n\n\n");


  logMessage( std::string("Catalog read in complete") );


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

    logMessage( std::string("Read 0 particles. Aborting") );

    printf(" Error reading particle file:%s\n", (userInput.getInputPart()).c_str() );
    printf("Specify the file using the variable \"partFile\", using the variable \"snapNum\", and making sure you are in the proper directory\n\n");
    exit(1);

  }


  // Saves this number of particles
  userInput.setNumParticles( numParticles );
  printf("\n Number of particles: %lli\n", userInput.getNumParticles() );


  // Allocate the partcles, only contains position information
  particlePosition *particle = new particlePosition[numParticles];
  printf(" Allocated %lli particles\n  Reading particle file...\n\n", numParticles );



  logMessage( std::string(    "Allocated " ) +
              std::to_string( numParticles ) +
              std::string(    " Particles" ) );



  // Read the particle positions
       numParticles  = readParticle( userInput, particle );

  if ( numParticles != userInput.getNumParticles() ){
    printf("Error: mismatch in particle read in\n Allocated: %lli\n Read in  : %lli\n", userInput.getNumParticles(), numParticles );

    logMessage( std::string("First  read in: ") + std::to_string( userInput.getNumParticles() ) + std::string(" particles") );
    logMessage( std::string("Second read in: ") + std::to_string(              numParticles   ) + std::string(" particles") );
    logMessage( std::string("Aborting."));

    exit(1);
  }


  logMessage( std::string(          "Read ") +
              std::to_string( numParticles ) +
              std::string(     " particles") );



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

  logMessage( std::string(                "Number of x cells= ") +
              std::to_string( (long long) userInput.getNrx()+1 ) );
  logMessage( std::string(                "Number of y cells= ") +
              std::to_string( (long long) userInput.getNry()+1 ) );
  logMessage( std::string(                "Number of x cells= ") +
              std::to_string( (long long) userInput.getNrz()+1 ) );
  logMessage( std::string(            "Total number of cells= ") +
              std::to_string( (long long) userInput.getNtotCell()  ) );


  std::cout << " Generating link list..." << std::endl;


  {
    logMessage( std::string( "Halo  memory: " ) + std::to_string( (long double) (sizeof( haloInfo         ) * userInput.getNumHalos()     / (1e6)) ) + std::string(" Mb") );
    logMessage( std::string( "Part  memory: " ) + std::to_string( (long double) (sizeof( particlePosition ) * userInput.getNumParticles() / (1e9)) ) + std::string(" Gb") );
    logMessage( std::string( "Link  memory: " ) + std::to_string( (long double) (sizeof( long long        ) * userInput.getNumParticles() / (1e9)) ) + std::string(" Gb") );
    logMessage( std::string( "Label memory: " ) + std::to_string( (long double) (sizeof( long long        ) * userInput.getNtotCell    () / (1e6)) ) + std::string(" Mb") );
    logMessage( std::string( "Total memory: " ) + std::to_string( (long double)((sizeof( long long        ) * userInput.getNtotCell    () +
                                                                                 sizeof( long long        ) * userInput.getNumParticles() +
                                                                                 sizeof( particlePosition ) * userInput.getNumParticles() +
                                                                                 sizeof( haloInfo         ) * userInput.getNumHalos()   ) / (1e9)) ) + std::string(" Gb") );
  }

  // Make link list between particles
  long long  *linkList = new long long [ userInput.getNumParticles() ];
  long long *labelList = new long long [ userInput.getNtotCell()     ];

  printf("  Allocated link  list of %lli elements\n", userInput.getNumParticles() );
  printf("  Allocated label list of %i elements\n"  , userInput.getNtotCell()     );


  logMessage( std::string(    "Allocated link list of "   ) +
              std::to_string( userInput.getNumParticles() ) +
              std::string(    " elements"                 ) );

  logMessage( std::string(                "Allocated label list of " ) +
              std::to_string(  (long long) userInput.getNtotCell()   ) +
              std::string(                 " elements"               ) );


  makeLinkList( userInput, particle, linkList, labelList );
  std::cout << " Done."     << std::endl  << std::endl;


  logMessage( std::string("Generated link list") );


  std::cout << " Refining halo list to those within our cells..." << std::endl;
  // Find the halos who's FOV is entirely inside the box
  {
    // Dummy array, useless
    haloInfo duhalos[2];

    N_halos = findBoxHalos( userInput, myHalos, duhalos, 0 ); //Returns the number of halos in the box on the first run through

    logMessage( std::string(    "Number of halos located in our box: " ) +
                std::to_string( (long long) N_halos                    ) );
  }


  haloInfo halos[N_halos];                                  //  Allocates our new halo array
  N_halos = findBoxHalos( userInput, myHalos,   halos, 1 ); //   Fills in our new halo array
  delete [] myHalos;                                        //        Deletes old halo array
  std::cout << "  Allocated " << N_halos    << " halos" << std::endl;


  logMessage( std::string(    "Allocated and filled new array of halos containing " ) +
              std::to_string( (long long) N_halos                                   ) +
              std::string(    " halos"                                              ) );


  std::cout << " Done."       << std::endl  << std::endl;


  logMessage( std::string( "New Halo  memory: " ) + std::to_string(  (long double) (sizeof( haloInfo ) * userInput.getNumHalos() / (1e6)) ) + std::string(" Mb") );


  /////////////////////////////////////////
  ////////Link halos, write files//////////
  /////////////////////////////////////////

  // Goes through link list, checking particles against halos
  // Locates particles in halo's Rvir
  //   box FOV x FOV x FOV
  //   box FOV x FOV x integration lengths
  // Writes image to FITS file

  std::cout << " Generating FITS images..." << std::endl << std::endl;
  linkHaloParticles( userInput, halos, particle, labelList, linkList );
  std::cout << " Done."                     << std::endl << std::endl;

  std::cout << "Wrote " << userInput.getNumFiles() << " Files " << std::endl;

  logMessage( std::string(       "Total number of output images: " ) +
              std::to_string( (long long) userInput.getNumFiles()  ) );

  delete[]  linkList;
  delete[] labelList;

  // Final log stuff
  // Record the run time in log


  int execution_end = clock();
  int runTime       = (execution_end - execution_start) / CLOCKS_PER_SEC;

  int         hours = runTime / 3600;
            runTime = runTime % 3600;
  int       minutes = runTime / 60;
            runTime = runTime % 60;

  logMessage( std::string(    "Elapsed runtime= " ) +
              std::to_string( (long long) hours   ) +
              std::string(    " hours, "          ) +
              std::to_string( (long long) minutes ) +
              std::string(    " minutes, "        ) +
              std::to_string( (long long) runTime ) +
              std::string(    " seconds"          ) );

}
