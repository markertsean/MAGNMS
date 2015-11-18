#include <iostream>
#include <fstream>
#include <cstring>
#include <sstream>
#include <iomanip>

#include <cmath>

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
      else if ( strcmp( inpC1, "useShortCatalog" ) == 0 ){
        myInput.setShortCat    ( std::stoi( inpS )  );
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
unsigned long readCatalog( haloInfo      halos[] ,  //Stores data of halos
                          inputInfo     userInfo ,  //Contains input data from user
                          unsigned long  N_halos ){ //0-just count valid halos, 1-store values



  //Generate short catalog name, based on input catalogs and flags
  int lastDot = (userInfo.getInputCatalog()).find_last_of(".");
  std::string shortCat = (userInfo.getInputCatalog()).substr( 0, lastDot + 1 );
  if ( userInfo.getShortCat() == 1 ){
    //Sets name based on flags
    if ( GLOBAL_MASS_LOWER_LIM != -1 ){
      std::stringstream stream;
      stream << std::fixed << std::setprecision(1) << log10( GLOBAL_MASS_LOWER_LIM );
      shortCat = shortCat + "Mmin_" + stream.str() +".";
    }
    if ( GLOBAL_MASS_UPPER_LIM != -1 ){
      std::stringstream stream;
      stream << std::fixed << std::setprecision(1) << log10( GLOBAL_MASS_UPPER_LIM );
      shortCat = shortCat + "Mmax_" + stream.str() +".";
    }
    shortCat = shortCat + "short_cat.DAT";
  }


  //First time for allocation, try to open file to see if it exists
  // if it exists, can read in first line to know size of halos
  // Otherwise, attempt to open normal file
  std::ifstream shortFile( shortCat );
  std::ifstream myFile( userInfo.getInputCatalog() );
  if ((              N_halos   == 0 ) &&  //First time read in
      ( userInfo.getShortCat() == 1 ) &&  //Using short catalog
      (       shortFile.good() == 1 ) ){  //Short file found
    userInfo.setCatType( "short" );
    std::cout << " Found short catalog: " << shortCat << std::endl;
  }
  else if( !myFile.is_open() ){
    std::cout << "\n Error opening halo catalog: " << userInfo.getInputCatalog() << std::endl;
    exit(1);
  }

  unsigned long N_read = 0;

  //Call read functions per function, to assign values, will use short if available
  if ( (userInfo.getCatType()).compare( "short" ) == 0 ){
    N_read = readShortCat ( shortFile, halos, N_halos );
  }
  else
  if ( (userInfo.getCatType()).compare(    "MD" ) == 0 ){
    N_read = readMultiDark(    myFile, halos, N_halos );
  }
  else
  {
      std::cout << " Unrecognized catalog type: " << userInfo.getCatType() << std::endl;
      exit(1);
  }

  //Write short catalog if one doesn't exist, and second time through
  if ( (userInfo.getShortCat() == 1) &&   //If using shortcat
               shortFile.bad() == 0  &&   //And couldn't find shortcat
       (             N_halos   >  0 ) ) { //And second time through
    std::ofstream writeFile( shortCat );
    if ( writeShortCat( writeFile, halos, N_halos ) ){
      std::cout << " Wrote short catalog: " << shortCat << std::endl;
    }
    else {
      std::cout << " Failed to write short catalog: " << shortCat << std::endl;
    }
  }

  return N_read;
}


//Reads short catalog for easy reading
unsigned long readShortCat ( std::ifstream &inpFile   ,
                             haloInfo         halos[] ,
                             unsigned long  N_halos   ){

  unsigned long num_Head;

  std::string str;  //The read in line

  //Reads header, which is one line with the number of halos
  // If we are just finding number to allocate, just read this
  std::getline( inpFile, str );
  std::stringstream num( str );
  num >> num_Head;
  if ( N_halos == 0 ){
    return num_Head;
  }

  float x, y, z, M, R, C, N, ba, ca;
  long id;
  long ds;

  int N_valid = 0;

  while ( std::getline( inpFile, str ) ){

    std::stringstream line( str );
    // x, y, z coordinate
    line >>    x;    line >>    y;     line >>    z;
    //M_tot, R_vir
    line >>    M;    line >>    R;
    // halo id number, Concentration, Number of halo, distinct/sub,
    line >>   id;    line >>    C;    line >>    N;    line >>   ds;

    // axis b/a ratio, c/a ratio, x axis, y axis, z axis
    line >>   ba;    line >>   ca;

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
        halos[ N_valid ].setBA( ba );
        halos[ N_valid ].setCA( ca );
        halos[ N_valid ].setDistinct( ds );
      }

      ++N_valid;
    }

  }


  return N_valid;
}


unsigned long readMultiDark( std::ifstream &inpFile   ,
                            haloInfo          halos[] ,
                            unsigned long   N_halos   ){

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
    line >>   ba;    line >>   ca;    line >> junk;    line >> junk;    line >> junk;

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
        halos[ N_valid ].setBA( ba );
        halos[ N_valid ].setCA( ca );
        halos[ N_valid ].setDistinct( ds );
      }

      ++N_valid;
    }

  }

  return N_valid;

}

//Writes short catalog for easy reading
bool         writeShortCat ( std::ofstream &inpFile   ,
                             haloInfo         halos[] ,
                             unsigned long  N_halos   ){

  //Writes number of valid halos at start

  inpFile << N_halos << std::endl;

  float x, y, z, M, R, C, N, ba, ca;
  long id;
  long ds;

  for ( int i=0; i<N_halos; ++i ){

    char writeLine[200];

    x  = halos[ i ].getX ();
    y  = halos[ i ].getY ();
    z  = halos[ i ].getZ ();
    C  = halos[ i ].getC ();
    M  = halos[ i ].getM ();
    N  = halos[ i ].getN ();
    R  = halos[ i ].getRm();
    id = halos[ i ].getID();
    ba = halos[ i ].getBA();
    ca = halos[ i ].getCA();
    ds = halos[ i ].getDistinct();

    sprintf(writeLine,"%12.4lf%12.4lf%12.4lf%12.4e%12.4lf%12lu%12.4lf%12.2lf%14lu%11.4lf%11.4lf\n",x,y,z,M,R,id,C,N,ds,ba,ca);

    inpFile << writeLine;
  }

  return true;
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
