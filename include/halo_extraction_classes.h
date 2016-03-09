#include <iostream>
#include <cstring>

#ifndef HE_CLASSES
#define HE_CLASSES

class inputInfo{

  public:

    inline inputInfo();
    static unsigned long Num_files;

    //Functions to change files the user may input
    void setReadFile    ( std::string   inpS ) { readFile        = inpS; }
    void setInputCatalog( std::string   inpS ) { inputCatalog    = inpS; }
    void setInputPart   ( std::string   inpS ) { inputPartFiles  = inpS; }
    void setInputHead   ( std::string   inpS ) { inputHeadFiles  = inpS; }
    void setParticleDir ( std::string   inpS ) { particleDir     = inpS; }
    void setHeaderDir   ( std::string   inpS ) {   headerDir     = inpS; }
    void setCatType     ( std::string   inpS ) {   catType       = inpS;
                                                  setPartMass();         }
    void setMinMass     ( double        inpD ) {   minMass       = inpD; }
    void setMaxMass     ( double        inpD ) {   maxMass       = inpD; }
    void setRadiusMult  ( double        inpD ) {radiusMultiplier = inpD; }
    void setShortCat    ( short         inpI ) {  useShortCat    = inpI; }
    void setUseShort    ( short         inpI ) {  usingShort     = inpI; }
    void setSnapNum     ( short         inpI ) {  snapshotNum    = inpI; }
    void setNumHalos    ( unsigned long inpI ) {  numHalos       = inpI; }
    void setNumParticles( long     long inpI ) {  numParticles   = inpI; }
    void setFOV         ( float         inpF ) {   boxFOV        = inpF; }
    void setXmin        ( float         inpF ) {   x_min         = inpF; }
    void setXmax        ( float         inpF ) {   x_max         = inpF; }
    void setYmin        ( float         inpF ) {   y_min         = inpF; }
    void setYmax        ( float         inpF ) {   y_max         = inpF; }
    void setZmin        ( float         inpF ) {   z_min         = inpF; }
    void setZmax        ( float         inpF ) {   z_max         = inpF; }
    void setNPixelsH    ( int           inpI ) {   N_pixels_h    = inpI; }
    void setNPixelsV    ( int           inpI ) {   N_pixels_v    = inpI; }
    void setIntegLength ( double        inpF ) {  maxIntegLength = inpF; }
    void setIntegStep   ( double        inpF ) {     integStep   = inpF; }



    std::string   getReadFile      () { return readFile         ; }
    std::string   getInputCatalog  () { return inputCatalog     ; }
    std::string   getInputPart     () { return inputPartFiles   ; }
    std::string   getInputHead     () { return inputHeadFiles   ; }
    std::string   getParticleDir   () { return particleDir      ; }
    std::string   getHeaderDir     () { return   headerDir      ; }
    std::string   getPartFileStart () { return partFileStart    ; }
    std::string   getCatType       () { return   catType        ; }

    double        getParticleMass  () { return   particleMass   ; }
    double        getMinMass       () { return        minMass   ; }
    double        getMaxMass       () { return        maxMass   ; }
    double        getRadiusMult    () { return radiusMultiplier ; }
    double        getFOV           () { return   boxFOV         ; }
    short         getShortCat      () { return   useShortCat    ; }
    short         getUsingShort    () { return   usingShort     ; }
    short         getSnapNum       () { return   snapshotNum    ; }

    unsigned long getNumHalos      () { return   numHalos       ; }
    long     long getNumParticles  () { return   numParticles   ; }

    float         getXmin          () { return   x_min          ; }
    float         getXmax          () { return   x_max          ; }
    float         getYmin          () { return   y_min          ; }
    float         getYmax          () { return   y_max          ; }
    float         getZmin          () { return   z_min          ; }
    float         getZmax          () { return   z_max          ; }
    float         getCell          () { return   cell           ; }

    int           getNlx           () { return   Nlx            ; }
    int           getNly           () { return   Nly            ; }
    int           getNlz           () { return   Nlz            ; }
    int           getNrx           () { return   Nrx            ; }
    int           getNry           () { return   Nry            ; }
    int           getNrz           () { return   Nrz            ; }
    int           getNtotCell      () { return   NtotCell       ; }

    int           getNPixlesH      () { return   N_pixels_h     ; }
    int           getNPixelsV      () { return   N_pixels_v     ; }

    double        getMaxIntegLength() { return   maxIntegLength ; }
    double        getIntegStep     () { return      integStep   ; }

    unsigned long getNumFiles      () { return        Num_files ; }



    // Keep track of number of times we wrote a file
    void wroteFile(){
      ++Num_files;
    }

    void setDirectory(){

      //Finds the last occurance of '/', which indicates directory
      int lastSlash = inputCatalog.find_last_of("/\\");
      catDir        = inputCatalog.substr( 0, lastSlash + 1 );
      catName       = inputCatalog.substr(    lastSlash + 1 );

      //Finds the last occurance of '.', which gives first file name
      int lastDot = catName.find_last_of(".");

      if ( catType.compare( "BMD") == 0 ){
        lastDot = ( catName.substr( 0, lastDot ) ).find_last_of(".");
      }

      //If directories for header and particle files
      // haven't been set, use catalog directory
      if ( particleDir=="" ){
          particleDir = catDir + catName.substr( 0, lastDot ) + "_MassMaps/";
      }

      if (   headerDir=="" ){
          headerDir   = catDir + catName.substr( 0, lastDot ) + "_Headers/";
      }

      partFileStart = catDir + "Part." ;
    }

    bool setBox( float inpCell ) {

      //If particle range hasn't been set, can't set the boxes
      if ( x_max == 0 || y_max == 0 || z_max == 0 ){
        return false;
      }

      cell = inpCell;

      Nlx = 0;
      Nly = 0;
      Nlz = 0;

      Nrx = ( x_max - x_min ) / cell ;     //Maximum box numbers
      Nry = ( y_max - y_min ) / cell ;
      Nrz = ( z_max - z_min ) / cell ;


      NtotCell = Nrx * Nry * Nrz;

      return true;
    }


    void setPartMass( ){

      if ( catType.compare( "BMDP"  ) == 0 ){ particleMass = 2.4e10; } else
      if ( catType.compare( "BMD"   ) == 0 ){ particleMass = 2.1e10; } else
      if ( catType.compare( "MDP"   ) == 0 ){ particleMass = 1.5e9 ; } else
      if ( catType.compare( "MD"    ) == 0 ){ particleMass = 8.7e9 ; } else
      if ( catType.compare( "BP"    ) == 0 ){ particleMass = 1.5e8 ; } else
      if ( catType.compare( "B"     ) == 0 ){ particleMass = 1.3e8 ; } else
      if ( catType.compare( "short" ) == 0 ){ /* No change needed */ } else
      {
        std::cout<<"Unrecognized catalog: "<<catType<<std::endl;
        exit(1);
      }
    }


  private:

    //User end input file
    std::string readFile      ; // Default extractInfo.dat, the file to read in user input, can be changed via command line argument

    //Data in input file, as to what files to use
    std::string inputCatalog  ; // Editable, the catalog of halos to use
    std::string inputPartFiles; // Editable, particle file to use
    std::string inputHeadFiles; // Editable, set the header file to use

    std::string particleDir   ; // Editable, location of the particle file
    std::string   headerDir   ; // Editable, location of particle header file
    std::string     catType   ; // Editable, default "MD", simulation catalog type, B, BP, MD, MDP, BMD, BMDP

    std::string  partFileStart; // Determine the start of the particle file, for passing to the fortran file
    std::string     catDir    ; // Directory of the halo catalog
    std::string     catName   ; // Name of the halo catalog


    // Halo restrictions imposed by user
    double          minMass    ; // Editable, minimum mass halo to use
    double          maxMass    ; // Editable, maximum mass halo to use
    double    radiusMultiplier ; // Editable, multiplier to radius to put in radius box, EDITEDITEDITEDITEDIT
    double              boxFOV ; // Editable, field of view of the image box

    short          useShortCat ; // Editable, use a short catalog, will generate a new catalog easier to read in
    short          snapshotNum ; // Editable, PMSS snapshot number to use

    int             N_pixels_h ; // Editable, number of pixels in each direction
    int             N_pixels_v ;

    double      maxIntegLength ; // Editable, maximum integration length in z, in Mpch
    double           integStep ; // Editable, step size in dex to use


    // Stuff set by the program
    unsigned long numHalos     ;
    long     long numParticles ;
    short         usingShort   ;

    double        particleMass ;  // Mass of a single particle, set by catalog type

    float                x_min ;  // Min/Max particle locations, used to find range of particles
    float                x_max ;
    float                y_min ;
    float                y_max ;
    float                z_min ;
    float                z_max ;

    float                 cell ; // Cell size in Mpc

    int                    Nlx ; // Number of cell boxes in each direction
    int                    Nly ;
    int                    Nlz ;
    int                    Nrx ;
    int                    Nry ;
    int                    Nrz ;
    int               NtotCell ;

};


//Default values initialized on construction
inputInfo::inputInfo( ){
  readFile       = "extractInfo.dat";
  inputCatalog   = "";
  inputPartFiles = "";
  inputHeadFiles = "";
  partFileStart  = "";
  particleDir    = "";
    headerDir    = "";
      catDir     = "";
      catName    = "";
      catType    = "MD";
       minMass    =  1e15;
       maxMass    = -1.0 ;
 radiusMultiplier =  1.0 ;
      useShortCat =  1   ;
      snapshotNum =  0   ;
      usingShort  =  0   ;

      boxFOV      =   8.0;//All Mpc

     numHalos     = 0;
     numParticles = 0;

            x_min = 0;
            x_max = 0;
            y_min = 0;
            y_max = 0;
            z_min = 0;
            z_max = 0;

            cell  = -1;

            Nlx   = -1;
            Nly   = -1;
            Nlz   = -1;
            Nrx   = -1;
            Nry   = -1;
            Nrz   = -1;
         NtotCell = -1;

       N_pixels_h = 1024;
       N_pixels_v = 1024;

     particleMass = 0.0;

          setPartMass();

   maxIntegLength = 400.0;
      integStep   = 0.2;

}


class haloInfo {

  public:
//14 params

    void setX ( float inpF ) {     x = inpF; }
    void setY ( float inpF ) {     y = inpF; }
    void setZ ( float inpF ) {     z = inpF; }
    void setC ( float inpF ) {     C = inpF; }
    void setM ( float inpF ) {     M = inpF; }
    void setN ( float inpI ) {     N = inpI; }
    void setID( long  inpL ) {    id = inpL; }
    void setBA( float inpF ) { ba_rat= inpF; }
    void setCA( float inpF ) { ca_rat= inpF; }

    void setRm( float inpF ) { float conversion = 1e-3; //kpc to Mpc, radius in kpc but all other coordinates Mpc
                                          R_max = inpF * conversion; }

    void setDistinct( long inpL ) { distinct = inpL; }


    float getX () { return  x     ; }
    float getY () { return  y     ; }
    float getZ () { return  z     ; }
    float getC () { return  C     ; }
    float getM () { return  M     ; }
    float getN () { return  N     ; }
    float getRm() { return  R_max ; }
    long  getID() { return  id    ; }
    float getBA() { return  ba_rat; }
    float getCA() { return  ca_rat; }

    long  getDistinct() { return  distinct; }


  private:

    float x       ;//=  0.0;  // x, y, z center of halo
    float y       ;//=  0.0;
    float z       ;//=  0.0;
    float C       ;//= -1.0;  // Concentration
    float M       ;//= -1.0;  // Mass of halo
    float N       ;//= -1  ;  // N_particles
    float R_max   ;//= -1.0;  // Rvir of halo
    long  id      ;//= -1  ;  // id number of halo
    long distinct ;//= -1  ;  // Distinct or subhalo, 0 or -1 d, #-distinct halo is subhalo of

    float ba_rat  ;//= -1.0;  // Ratio of axis b to a
    float ca_rat  ;//= -1.0;  // Ratio of axis c to a

};


class particlePosition{

public:

    inline particlePosition();

    float x_pos, y_pos, z_pos;

};

particlePosition::particlePosition(){
  x_pos = 0;
  y_pos = 0;
  z_pos = 0;
}


extern "C" long long readpmss_( int *jstep, char *filestart, int *filestartlength );

#endif // HE_CLASSES
