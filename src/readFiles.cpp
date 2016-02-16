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
      else if ( strcmp( inpC1, "snapNum" ) == 0 ){
        myInput.setSnapNum     ( std::stoi( inpS )  );
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
             !shortFile.good()       &&   //And couldn't find shortcat
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



//Locates the particle file, or generates one from PMSS
long long setPartFile( inputInfo &userInput ) { //All the user input


  long long numParticles = 0;


  //If the user specified a file, use it
  if ( (userInput.getInputPart()!="") ){

    std::string partCatalog = userInput.getInputPart();

    char tempS[60];

    sprintf(tempS,"%s.header.dat", (partCatalog.substr(0,partCatalog.length()-4)).c_str() );


    //Attempts to open the files, to see if they exist
    std::ifstream testPart( partCatalog );
    std::ifstream testHead(       tempS );


    //If we could open both files, that's all
    if ( (testPart.good() == 1) && (testHead.good() == 1) ){

      printf( "%s%s\n", " Using particle file: ", ( userInput.getInputPart() ).c_str() );
      printf( "%s%s\n", " Using   header file: ",                                tempS );

      //Go through header file to determine number of particles for array allocation
      std::string   tempStr[8];
      double tempInt[8];

      for ( int i = 0; i <= 7; ++i ){ //Number of particles is the 8th item
        testHead >> tempStr[i];
        testHead >> tempInt[i];
      }
      numParticles = tempInt[7];

      userInput.setInputHead( tempS );

      return numParticles;
    }
  }

  //If nothing specified, or specified file not found

  int snapNum = userInput.getSnapNum();                     //Snapshot number
  std::string partFileStart = userInput.getPartFileStart(); //Beginning of particle file name


  //Generates file names for particle file, header file, based on beginning of file name
  char testPartFile[60];
  char testHeadFile[60];
  sprintf(testPartFile,"%s%4.4i.DAT"       , partFileStart.c_str(), snapNum );
  sprintf(testHeadFile,"%s%4.4i.header.dat", partFileStart.c_str(), snapNum );


  //Attempts to open the files, to see if they exist
  std::ifstream testPart( testPartFile );
  std::ifstream testHead( testHeadFile );


  //If we could open both files, that's all
  if ( (testPart.good() == 1) && (testHead.good() == 1) ){

    printf( "%s%s\n", " Using particle file: ", testPartFile );
    printf( "%s%s\n", " Using   header file: ", testHeadFile );

    //Go through header file to determine number of particles for array allocation
    std::string   tempStr[8];
    double tempInt[8];

    for ( int i = 0; i <= 7; ++i ){ //Number of particles is the 8th item
      testHead >> tempStr[i];
      testHead >> tempInt[i];
    }
    numParticles = tempInt[7];

    userInput.setInputPart( testPartFile );
    userInput.setInputHead( testHeadFile );

    return numParticles;

  }
  //If we couldn't open both files, read in from pmss file using other input data
  else{

    char tempC[60];                                           //For passing begining of file name to the fortran function
                                                                // if blank, just uses current file directory
    strcpy(tempC, partFileStart.c_str());                     //Copy the string to the char array
    int fileNameLength = partFileStart.length();              //Length of the particle file name, needed for
                                                                // fortran format string
    /*
      readpmss reads the fortran unformatted binary file for the particles
      file name from the snapshop
      and file location passed to it
      Outputs a header file with basic info,
      and a shitton of particles in a file
      now need to read in from the C end,
      allocate the location arrays,
      and then scan used halo catalog for the acceptance criteria
      Also, it returns the number of particles in the file
    */
    numParticles = readpmss_( &snapNum, tempC, &fileNameLength );

    userInput.setInputPart( testPartFile );
    userInput.setInputHead( testHeadFile );

    return numParticles;

  }

  return 0;
}



bool readParticle( inputInfo userInfo, particlePosition particles[] ){



  std::ifstream inputParticleFile( ( userInfo.getInputPart() ).c_str() );

  if ( inputParticleFile.good() ){
    for ( int i = 0; i < userInfo.getNumParticles() ; ++i ){

      inputParticleFile >> particles[i].x_pos;
      inputParticleFile >> particles[i].y_pos;
      inputParticleFile >> particles[i].z_pos;

    }
  }
  else {

    printf(" Error opening particle file: %s\n\n", ( userInfo.getInputPart() ).c_str() );
    exit(1);

  }

  return true;
}
