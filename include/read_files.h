#include <cstring>

#ifndef READ_FILES
#define READ_FILES

  //Reads input files and optional directories to write files to (otherwise default)
  bool readUserInput( std::string fileName, inputInfo &myInput );

  //Reads halo catalog, and saves in halos array
  unsigned long readCatalog( haloInfo   halos[] ,  //Stores data of halos
                             inputInfo userInfo ,  //Contains input data from user
                             int        N_halos ); //0-just count valid halos, 1-store values

  //Reads catalog with Multidark format
  unsigned long readMultiDark( std::ifstream &inpFile   ,
                               haloInfo         halos[] ,
                               int            N_halos   );

  //Test halos against user constraints
  bool validHalo( float  inpM ,  //Input halos mass
                  long  inpDS ); //Distinct/sub flag

#endif
