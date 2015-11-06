#include <iostream>
#include <cstring>

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


  //Stores the user input from the input file
  inputInfo userInput;


  //User can specify input file via command line, default: extractInfo.dat
  if ( arg > 1 ){
    userInput.setReadFile( argv[1] );
  }
  printf("\nUsing input file: %s\n", (userInput.getReadFile()).c_str());


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


  //Have input files and directories to write to.
  //Need to:
  //        Include halo selection criteria, mass, distinct, coordinates,
  //                radius, triaxiality, ...
  //        Create directories
  //        Write overarching header file, with used halos
  //        Read in short halo list, if catalog already used
  //        Read in particles, maybe use link list to speed up process?




}
