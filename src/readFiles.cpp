#include <iostream>
#include <fstream>
#include <cstring>
#include <sstream>
#include "halo_extraction_classes.h"
#include "read_files.h"


double GLOBAL_MASS_UPPER_LIM = -1;
double GLOBAL_MASS_LOWER_LIM = -1;


//Reads input files and optional directories to write files to (otherwise default)
bool readUserInput( std::string fileName, inputInfo &myInput ){

  FILE *pFile;
  char   inpC1[50],inpC2[50];
  pFile = fopen( fileName.c_str() ,"r");

  //If we can open the file, cycle through lines
  if ( pFile != NULL ){
    while ( fscanf(pFile,"%s%s",inpC1,inpC2) != EOF ){

      std::string inpS = std::string(inpC2);

      if      ( strcmp( inpC1, "haloFile" ) == 0 ){
        myInput.setInputCatalog( inpS  );
      }
      else if ( strcmp( inpC1, "partFile" ) == 0 ){
        myInput.setInputPart   ( inpS  );
      }
      else if ( strcmp( inpC1, "headDir"  ) == 0 ){
        myInput.setParticleDir ( inpS  );
      }
      else if ( strcmp( inpC1, "partDir"  ) == 0 ){
        myInput.setHeaderDir   ( inpS  );
      }
      else if ( strcmp( inpC1, "catType"  ) == 0 ){
        myInput.setHeaderDir   ( inpS  );
      }
      else if ( strcmp( inpC1, "minMass"  ) == 0 ){
        myInput.setMinMass     ( std::stod( inpS )  );
      }
      else if ( strcmp( inpC1, "maxMass"  ) == 0 ){
        myInput.setMaxMass     ( std::stod( inpS )  );
      }
      else if ( strcmp( inpC1, "rMult"    ) == 0 ){
        myInput.setRadiusMult  ( std::stod( inpS )  );
      }
    }
  }
  //If file doesn't open, throw error
  else {
    std::cout << "Could not open input file: " << fileName << std::endl;
    return false;
  }

  fclose(pFile);

  if ( myInput.getInputCatalog() == "" || myInput.getInputPart() == "" ){
    printf("Missing required input files:\n   Halo catalog (haloFile) = %s\n   Particle file (partFile) = %s\n\n",
          (myInput.getInputCatalog()).c_str(),(myInput.getInputPart()).c_str());
    return false;
  }

  myInput.setDirectory();

  GLOBAL_MASS_LOWER_LIM = myInput.getMinMass();
  GLOBAL_MASS_UPPER_LIM = myInput.getMaxMass();

  return true;
}


//Reads halo catalog, and saves in halos array
unsigned long readCatalog( haloInfo   halos[] ,  //Stores data of halos
                          inputInfo userInfo ,  //Contains input data from user
                          int        N_halos ){ //0-just count valid halos, 1-store values


  std::ifstream myFile( userInfo.getInputCatalog() );
  if( !myFile.is_open() ){
    std::cout << "\n Error opening halo catalog: " << userInfo.getInputCatalog() << std::endl;
    exit(1);
  }

  //Call read functions per function, to assign values
  //Different function if it's a file we used, will write to new file
  N_halos = readMultiDark( myFile, halos, N_halos );

  return N_halos;
}


unsigned long readMultiDark( std::ifstream &inpFile   ,
                            haloInfo         halos[] ,
                            int            N_halos   ){

  //Skip the header of file
  int headerLength = 17;
  {
    std::string junk;
    for (int i=0;i<headerLength;++i)
      std::getline(inpFile, junk);
  }


  std::string str;  //The read in line
  std::string junk; //For stuff we don't need

  float   x, y, z, M, R, C, N, ba, ca, xa, ya, za;
  long id;
  long ds;

  int N_valid = 0;

  while ( std::getline( inpFile, str ) ){

    std::stringstream line( str );
    // x, y, z coordinate
    line >>    x;    line >>    y;     line >>    z;    line >> junk;    line >> junk;     line >> junk;    line >> junk;

    //M_tot, R_vir
    line >>    M;    line >>    R;    line >> junk;    line >> junk;

    // halo id number, Concentration, Number of halo, distinct/sub,
    line >>   id;    line >>    C;    line >>    N;    line >>   ds;    line >> junk;    line >> junk;    line >> junk;    line >> junk;

    // axis b/a ratio, c/a ratio, x axis, y axis, z axis
    line >>   ba;    line >>   ca;    line >>   xa;    line >>   ya;    line >>   za;

    //Test if halo is distinct and in mass range
    if ( validHalo ( M, ds ) ){

      //Can only occur on second run
      //First run returns number of valid halos
      //halos array then allocated
      //Then run through and save those values
      if ( N_halos > 0 ){
        halos[ N_valid ].setX (  x );
        halos[ N_valid ].setY (  y );
        halos[ N_valid ].setZ (  z );
        halos[ N_valid ].setC (  C );
        halos[ N_valid ].setM (  M );
        halos[ N_valid ].setN (  N );
        halos[ N_valid ].setRm(  R );
        halos[ N_valid ].setID( id );
        halos[ N_valid ].setXa( xa );
        halos[ N_valid ].setYa( ya );
        halos[ N_valid ].setZa( za );
        halos[ N_valid ].setBA( ba );
        halos[ N_valid ].setCA( ca );
        halos[ N_valid ].setDistinct( ds );
      }

      ++N_valid;
    }

  }

  return N_valid;

}


//Tests for validity of halo, based on halo mass and distinct vs subhalo
bool validHalo( float inpM, long inpDS ){

  //Using an upper limit, and mass of halo above the upper limit
  if ( (GLOBAL_MASS_UPPER_LIM > -1) && ( inpM > GLOBAL_MASS_UPPER_LIM ) )
    return false;

  //Using an lower limit, and mass of halo below the lower limit
  if ( (GLOBAL_MASS_LOWER_LIM > -1) && ( inpM < GLOBAL_MASS_LOWER_LIM ) )
    return false;

  //If a subhalo, toss it. DS >0
  if ( inpDS >0 )
    return false;

  return true;
}
