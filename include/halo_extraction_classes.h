#include <iostream>
#include <cstring>

#ifndef HE_CLASSES
#define HE_CLASSES

class inputInfo{

  public:

    //Functions to change files the user may input
    void setReadFile    ( std::string inpS ) { readFile       = inpS; }
    void setInputCatalog( std::string inpS ) { inputCatalog   = inpS; }
    void setInputPart   ( std::string inpS ) { inputPartFiles = inpS; }
    void setParticleDir ( std::string inpS ) { particleDir    = inpS; }
    void setHeaderDir   ( std::string inpS ) {   headerDir    = inpS; }

    std::string getReadFile    () { return readFile       ; }
    std::string getInputCatalog() { return inputCatalog   ; }
    std::string getInputPart   () { return inputPartFiles ; }
    std::string getParticleDir () { return particleDir    ; }
    std::string getHeaderDir   () { return   headerDir    ; }

    void setDirectory(){

      //Finds the last occurance of '/', which indicates directory
      int lastSlash = inputCatalog.find_last_of("/\\");
      catDir        = inputCatalog.substr( 0, lastSlash + 1 );
      catName       = inputCatalog.substr(    lastSlash + 1 );

      //Finds the last occurance of '.', which gives first file name
      int lastDot = catName.find_last_of(".");

      //If directories for header and particle files
      // haven't been set, use catalog directory
      if ( particleDir=="" ){
        particleDir = catDir + catName.substr( 0, lastDot ) + "_Particles/";
      }

      if (   headerDir=="" ){
          headerDir = catDir + catName.substr( 0, lastDot ) + "_Headers/";
      }

    }

  private:
    //User end input file
    std::string readFile       = "extractInfo.dat";

    //Data in input file, as to what files to use
    std::string inputCatalog   = "";
    std::string inputPartFiles = "";

    std::string particleDir    = "";
    std::string   headerDir    = "";

    std::string     catDir     = "";
    std::string     catName    = "";

    //Halo restrictions imposed by user
    double          minMass    =  1e15;
    double          maxMass    = -1.0 ;
    double    radiusMultiplier =  1.0 ;
};


class haloInfo {




  private:
    float xp[3] = {0,0,0};
    float M     = -1.0;
    float R_max = -1.0;
    float C     = -1.0;
    long  id    =    0;

};

#endif // HE_CLASSES
