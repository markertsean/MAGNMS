#include <iostream>
#include <cstring>

#ifndef HE_CLASSES
#define HE_CLASSES

class inputInfo{

  public:

    inline inputInfo();

    //Functions to change files the user may input
    void setReadFile    ( std::string inpS ) { readFile        = inpS; }
    void setInputCatalog( std::string inpS ) { inputCatalog    = inpS; }
    void setInputPart   ( std::string inpS ) { inputPartFiles  = inpS; }
    void setParticleDir ( std::string inpS ) { particleDir     = inpS; }
    void setHeaderDir   ( std::string inpS ) {   headerDir     = inpS; }
    void setCatType     ( std::string inpS ) {   catType       = inpS; }
    void setMinMass     ( double      inpD ) {   minMass       = inpD; }
    void setMaxMass     ( double      inpD ) {   maxMass       = inpD; }
    void setRadiusMult  ( double      inpD ) {radiusMultiplier = inpD; }
    void setShortCat    ( short       inpI ) {  useShortCat    = inpI; }

    std::string getReadFile    () { return readFile         ; }
    std::string getInputCatalog() { return inputCatalog     ; }
    std::string getInputPart   () { return inputPartFiles   ; }
    std::string getParticleDir () { return particleDir      ; }
    std::string getHeaderDir   () { return   headerDir      ; }
    std::string getCatType     () { return   catType        ; }
    double      getMinMass     () { return   minMass        ; }
    double      getMaxMass     () { return   maxMass        ; }
    double      getRadiusMult  () { return radiusMultiplier ; }
    short       getShortCat    () { return   useShortCat    ; }

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
    std::string readFile      ;// = "extractInfo.dat";

    //Data in input file, as to what files to use
    std::string inputCatalog  ;// = "";
    std::string inputPartFiles;// = "";

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
};

//Default values initialized on construction
inputInfo::inputInfo( ){
  readFile       = "extractInfo.dat";
  inputCatalog   = "";
  inputPartFiles = "";
  particleDir    = "";
    headerDir    = "";
      catDir     = "";
      catName    = "";
      catType    = "MD";
       minMass    =  1e15;
       maxMass    = -1.0 ;
 radiusMultiplier =  1.0 ;
      useShortCat =  1   ;

      boxFOV      =   8.0;//All Mpc
      losLength1  =  50.0;
      losLength2  = 250.0;
      losLength3  = 500.0;
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
    void setRm( float inpF ) { R_max = inpF; }
    void setID( long  inpL ) {    id = inpL; }
    void setBA( float inpF ) { ba_rat= inpF; }
    void setCA( float inpF ) { ca_rat= inpF; }
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


/*
extern "C" {
               extern struct{
                   double x;
                   int a, b, c;
               } abc_;
             }
*/

/*
Integer*4  ::    nx, ny, nz, Nbuffer                               ! number of domains
Real*4     ::     xL,xR,yL,yR,zL,zR,dBuffer   ! boundaries of domain
                     ! ------------------------- cosmology ----------------------
Real*4     ::  AEXPN, Om0,Oml0,hubble        ! cosmology
Real*4     ::  Box,MassOne,fraction                   !  current simulation
Integer*4  ::  jStep
Integer*4  ::  Np,Ntot                          ! Number of particles

                     ! -------------------------  Particles ----------------------

Real*4,       ALLOCATABLE,   DIMENSION(:) ::             &
                                      Xpp,  Ypp,  Zpp!,   &      !   coords
!                                      Wxp,  Wyp,  Wzp
Integer*8       ::  iCount
*/

extern "C" {
extern struct{
int nx, ny, nz;
float xl, xr, yl, yr, zl, zr, dBuffer;
float aexpn, om0, oml0, hubble;
float box,massone,fraction;
int jstep;
int np,ntot;
float* xpp;
float* ypp;
float* zpp;
int icount;
void readconfig_();
} setarrs_;
}

extern "C" void teststuff_();

#endif // HE_CLASSES
