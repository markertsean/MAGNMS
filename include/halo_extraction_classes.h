#include <iostream>
#include <cstring>

#ifndef HE_CLASSES
#define HE_CLASSES

class inputInfo{

  public:

    inline inputInfo();

    //Functions to change files the user may input
    void setReadFile    ( std::string   inpS ) { readFile        = inpS; }
    void setInputCatalog( std::string   inpS ) { inputCatalog    = inpS; }
    void setInputPart   ( std::string   inpS ) { inputPartFiles  = inpS; }
    void setInputHead   ( std::string   inpS ) { inputHeadFiles  = inpS; }
    void setParticleDir ( std::string   inpS ) { particleDir     = inpS; }
    void setHeaderDir   ( std::string   inpS ) {   headerDir     = inpS; }
    void setCatType     ( std::string   inpS ) {   catType       = inpS; }
    void setMinMass     ( double        inpD ) {   minMass       = inpD; }
    void setMaxMass     ( double        inpD ) {   maxMass       = inpD; }
    void setRadiusMult  ( double        inpD ) {radiusMultiplier = inpD; }
    void setShortCat    ( short         inpI ) {  useShortCat    = inpI; }
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



    std::string   getReadFile      () { return readFile         ; }
    std::string   getInputCatalog  () { return inputCatalog     ; }
    std::string   getInputPart     () { return inputPartFiles   ; }
    std::string   getInputHead     () { return inputHeadFiles   ; }
    std::string   getParticleDir   () { return particleDir      ; }
    std::string   getHeaderDir     () { return   headerDir      ; }
    std::string   getPartFileStart () { return partFileStart    ; }
    std::string   getCatType       () { return   catType        ; }

    double        getMinMass       () { return   minMass        ; }
    double        getMaxMass       () { return   maxMass        ; }
    double        getRadiusMult    () { return radiusMultiplier ; }
    double        getFOV           () { return   boxFOV         ; }
    short         getShortCat      () { return   useShortCat    ; }
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

      partFileStart = catDir + "Part." ;//+ std::to_string( snapshotNum ) + ""

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

  private:
    //User end input file
    std::string readFile      ;// = "extractInfo.dat";

    //Data in input file, as to what files to use
    std::string inputCatalog  ;// = "";
    std::string inputPartFiles;// = "";
    std::string inputHeadFiles;

    std::string  partFileStart;

    std::string particleDir   ;// = "";
    std::string   headerDir   ;// = "";

    std::string     catDir    ;// = "";
    std::string     catName   ;// = "";

    std::string     catType   ;// = "MD";

    //Halo restrictions imposed by user
    double          minMass    ;//=  1e15;
    double          maxMass    ;//= -1.0 ;
    double    radiusMultiplier ;//=  1.0 ;

    double              boxFOV ;
    double          losLength1 ;
    double          losLength2 ;
    double          losLength3 ;

    short          useShortCat ;//=  1   ;

    short          snapshotNum ;

    unsigned long numHalos     ;
    long     long numParticles ;

    float x_min, x_max, y_min, y_max, z_min, z_max;

    float cell;

    int Nlx, Nly, Nlz, Nrx, Nry, Nrz, NtotCell;

    int N_pixels_h, N_pixels_v;

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

      boxFOV      =   8.0;//All Mpc
      losLength1  =  50.0;
      losLength2  = 250.0;
      losLength3  = 500.0;

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


extern "C" unsigned long readpmss_( int *jstep, char *filestart, int *filestartlength );

#endif // HE_CLASSES
