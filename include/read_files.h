#include <cstring>


#ifndef READ_FILES
#define READ_FILES

  //Reads input files and optional directories to write files to (otherwise default)
  bool readUserInput( std::string fileName, inputInfo &myInput );

  //Reads halo catalog, and saves in halos array
  unsigned long readCatalog( haloInfo      halos[] ,  //Stores data of halos
                             inputInfo   *userInfo ,  //Contains input data from user
                             unsigned long N_halos ); //0-just count valid halos, 1-store values

  //Reads short catalog for easy reading
  unsigned long readShortCat ( std::ifstream &inpFile   ,
                               haloInfo         halos[] ,
                               unsigned long  N_halos   );


  //Reads catalog with Multidark format
  unsigned long readMultiDark( std::ifstream &inpFile   ,
                               haloInfo         halos[] ,
                               unsigned long  N_halos   );

  //Reads catalog with Multidark format
  unsigned long readBigMultiDark( std::string    inpFile   ,
                                  haloInfo         halos[] ,
                                  unsigned long  N_halos   );

  //Writes short catalog for easy reading
  bool         writeShortCat ( std::ofstream &inpFile   ,
                               haloInfo         halos[] ,
                               unsigned long  N_halos   );

  //Test halos against user constraints
  bool validHalo( float  inpM ,  //Input halos mass
                  long  inpDS ); //Distinct/sub flag

  //Locates the particle file, or generates one from PMSS
  long long setPartFile( inputInfo &userInfo ) ; //All the user input

  //Reads the particle file
  long long readParticle( inputInfo userInfo, particlePosition *particles );

  //Reads the header associated with the file
  bool readHeader  ( inputInfo &userInfo );


#endif
