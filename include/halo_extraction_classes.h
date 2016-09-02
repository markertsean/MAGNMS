#include <iostream>
#include <cstring>
#include <fstream>

#ifndef HE_CLASSES
#define HE_CLASSES



// Log file name
extern std::string logFileName;

// Generates log files, printing text to file
inline void logMessage( const std::string &text ){

    std::ofstream log_file(  logFileName, std::ios_base::out | std::ios_base::app );
    log_file << text << std::endl;

}



class inputInfo{

  public:

    inline inputInfo();
    static unsigned long Num_files;

    //Functions to change files the user may input
    void setReadFile     ( std::string   inpS ) { readFile        = inpS; }
    void setInputCatalog ( std::string   inpS ) { inputCatalog    = inpS; }
    void setInputPart    ( std::string   inpS ) { inputPartFiles  = inpS; }
    void setInputHead    ( std::string   inpS ) { inputHeadFiles  = inpS; }
    void setParticleDir  ( std::string   inpS ) { particleDir     = inpS; }
    void setHeaderDir    ( std::string   inpS ) {   headerDir     = inpS; }
    void setOutputDir    ( std::string   inpS ) {   outputDir     = inpS; }
    void setCatType      ( std::string   inpS ) {   catType       = inpS;
                                                    setPartMass();        }

    void setIntegAxis    ( char          inpC ) {  if ( inpC != 'x' &&
                                                        inpC != 'y' &&
                                                        inpC != 'z') {

                                                      std::cout << "Unrecognized axis: " << inpC << std::endl;
                                                      exit(1);
                                                    }
                                                   integAxis      = inpC; }

    void setMinMass      ( double        inpD ) {   minMass       = inpD; }
    void setMaxMass      ( double        inpD ) {   maxMass       = inpD; }
    void setRadiusMult   ( double        inpD ) {radiusMultiplier = inpD; }
    void setPhysicalSize ( double        inpD ) { physicalSize    = inpD; }
    void setRedshift     ( double        inpD ) { redshift        = inpD; }
    void setRadiusConvert( double        inpD ) {radiusConverter  = inpD; }
    void setShortCat     ( short         inpI ) {  useShortCat    = inpI; }
    void setUseShort     ( short         inpI ) {  usingShort     = inpI; }
    void setSnapNum      ( short         inpI ) {  snapshotNum    = inpI; }
    void setNumHalos     ( unsigned long inpI ) {  numHalos       = inpI; }
    void setNumParticles ( long     long inpI ) {  numParticles   = inpI; }
    void setFOV          ( float         inpF ) {   boxFOV        = inpF; }
    void setXmin         ( float         inpF ) {   x_min         = inpF; }
    void setXmax         ( float         inpF ) {   x_max         = inpF; }
    void setYmin         ( float         inpF ) {   y_min         = inpF; }
    void setYmax         ( float         inpF ) {   y_max         = inpF; }
    void setZmin         ( float         inpF ) {   z_min         = inpF; }
    void setZmax         ( float         inpF ) {   z_max         = inpF; }
    void setNPixelsH     ( int           inpI ) {   N_pixels_h    = inpI; }
    void setNPixelsV     ( int           inpI ) {   N_pixels_v    = inpI; }
    void setIntegLength  ( double        inpF ) {  maxIntegLength = inpF; }
    void setIntegStep    ( double        inpF ) {     integStep   = inpF; }
    void setPMssFirstNum ( int           inpI ) { PMssFirstFileNum= inpI; }
    void setPMssLastNum  ( int           inpI ) { PMssLastFileNum = inpI; }

    void setAexpn        ( float         inpF ) {  aexpn          = inpF; }
    void setHubble       ( float         inpF ) {  hubble         = inpF; }
    void setOmega_m      ( float         inpF ) {  om_o           = inpF; }
    void setOmega_l      ( float         inpF ) {  om_l           = inpF; }
    void setParticleMass ( double        inpF ) {  particleMass   = inpF; }


    std::string   getReadFile      () const { return readFile         ; }
    std::string   getInputCatalog  () const { return inputCatalog     ; }
    std::string   getInputPart     () const { return inputPartFiles   ; }
    std::string   getInputHead     () const { return inputHeadFiles   ; }
    std::string   getParticleDir   () const { return particleDir      ; }
    std::string   getHeaderDir     () const { return   headerDir      ; }
    std::string   getPartFileStart () const { return partFileStart    ; }
    std::string   getCatType       () const { return   catType        ; }
    std::string   getOutputDir     () const { return   outputDir      ; }

    char          getIntegAxis     () const { return   integAxis      ; }

    double        getParticleMass  () const { return   particleMass   ; }
    double        getMinMass       () const { return        minMass   ; }
    double        getMaxMass       () const { return        maxMass   ; }
    double        getRadiusMult    () const { return radiusMultiplier ; }
    double        getPhysicalSize  () const { return     physicalSize ; }
    double        getPixelUnit     () const { return     pixelUnits   ; }
    double        getRedshift      () const { return     redshift     ; }
    double        getRadiusConvert () const { return radiusConverter  ; }
    double        getFOV           () const { return   boxFOV         ; }
    short         getShortCat      () const { return   useShortCat    ; }
    short         getUsingShort    () const { return   usingShort     ; }
    short         getSnapNum       () const { return   snapshotNum    ; }


    unsigned long getNumHalos      () const { return   numHalos       ; }
    long     long getNumParticles  () const { return   numParticles   ; }

    float         getXmin          () const { return   x_min          ; }
    float         getXmax          () const { return   x_max          ; }
    float         getYmin          () const { return   y_min          ; }
    float         getYmax          () const { return   y_max          ; }
    float         getZmin          () const { return   z_min          ; }
    float         getZmax          () const { return   z_max          ; }
    float         getCell          () const { return   cell           ; }

    int           getNlx           () const { return   Nlx            ; }
    int           getNly           () const { return   Nly            ; }
    int           getNlz           () const { return   Nlz            ; }
    int           getNrx           () const { return   Nrx            ; }
    int           getNry           () const { return   Nry            ; }
    int           getNrz           () const { return   Nrz            ; }
    int           getNtotCell      () const { return   NtotCell       ; }

    int           getNPixlesH      () const { return   N_pixels_h     ; }
    int           getNPixelsV      () const { return   N_pixels_v     ; }

    int           getPMssFirstNum  () const { return PMssFirstFileNum ; }
    int           getPMssLastNum   () const { return PMssLastFileNum  ; }

    double        getMaxIntegLength() const { return   maxIntegLength ; }
    double        getIntegStep     () const { return      integStep   ; }

    unsigned long getNumFiles      () const { return        Num_files ; }

    float         getAexpn         () const { return   aexpn          ; }
    float         getHubble        () const { return   hubble         ; }
    float         getOmega_m       () const { return   om_o           ; }
    float         getOmega_l       () const { return   om_l           ; }


    // Keep track of number of times we wrote a file
    void wroteFile(){
      ++Num_files;
    }

    void setDirectory(){

      //Finds the last occurance of '/', which indicates directory
      int lastSlashCat =   inputCatalog.find_last_of("/\\");
      int lastSlashPar = inputPartFiles.find_last_of("/\\");
      catDir        = inputCatalog.substr( 0, lastSlashCat + 1 );
      catName       = inputCatalog.substr(    lastSlashCat + 1 );

      //Finds the last occurance of '.', which gives first file name
      int lastDot = catName.find_last_of(".");

      if ( catType.compare( "BMD") == 0 || catType.compare( "MDP") == 0 ){
        lastDot = ( catName.substr( 0, lastDot ) ).find_last_of(".");
      }

      //If directories for header and particle files
      // haven't been set, use catalog directory
//      if ( particleDir=="" ){
//          particleDir = particleDir + catName.substr( 0, lastDot ) + "_MassMaps/";
//      }

      if (   headerDir=="" ){
          headerDir   = catDir      + catName.substr( 0, lastDot ) + "_Headers/";
      }

      partFileStart = particleDir + "Part." ;
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

      Nrx = ( x_max - x_min ) / cell;     //Maximum box numbers
      Nry = ( y_max - y_min ) / cell;
      Nrz = ( z_max - z_min ) / cell;


      NtotCell = (Nrx+1) * (Nry+1) * (Nrz+1);

      return true;
    }

    float getCatBoxSize(){
      if ( catType.compare( "B"   ) == 0 || catType.compare(   "BP" ) == 0 ){
        return  250.0;
      }
      if ( catType.compare( "MD"  ) == 0 || catType.compare(  "MDP" ) == 0 ){
        return 1000.0;
      }
      // ( catType.compare( "BMD" ) == 0 || catType.compare( "BMDP" ) == 0 )
        return 2500.0;
    }


    void setPartMass( ){

      if ( catType.compare( "BMDP"  ) == 0 ){ pixelUnits = 0.6777; } else
      if ( catType.compare( "BMD"   ) == 0 ){ pixelUnits = 0.7100; } else
      if ( catType.compare( "MDP"   ) == 0 ){ pixelUnits = 0.6777; } else
      if ( catType.compare( "MD"    ) == 0 ){ pixelUnits = 0.7100; } else
      if ( catType.compare( "BP"    ) == 0 ){ pixelUnits = 0.6777; } else
      if ( catType.compare( "B"     ) == 0 ){ pixelUnits = 0.7100; } else
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
    std::string  outputDir    ;


    // Halo restrictions imposed by user
    double          minMass    ; // Editable, minimum mass halo to use
    double          maxMass    ; // Editable, maximum mass halo to use
    double    radiusMultiplier ; // Editable, multiplier to radius to put in radius box, EDITEDITEDITEDITEDIT
    double    radiusConverter  ; // Editable, convert units to match the position coordinates
    double              boxFOV ; // Editable, field of view of the image box

    short          useShortCat ; // Editable, use a short catalog, will generate a new catalog easier to read in
    short          snapshotNum ; // Editable, PMSS snapshot number to use

    int             N_pixels_h ; // Editable, number of pixels in each direction
    int             N_pixels_v ;

    double        physicalSize ; // In degrees, required for FITS header
    double          redshift   ; // Required for fits header

    double      maxIntegLength ; // Editable, maximum integration length in z, in Mpch
    double           integStep ; // Editable, step size in dex to use


    // Stuff set by the program
    unsigned long numHalos     ;
    long     long numParticles ;
    short         usingShort   ;

    double        particleMass ;  // Mass of a single particle, set by catalog type
    double        pixelUnits   ;

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

    // Info from the simulation

    float aexpn;
    float om_o;
    float om_l;
    float hubble;



    ///////////////////////////////////////
    ///////////////////////////////////////
    //////        NEW           ///////////
    ///////////////////////////////////////
    ///////////////////////////////////////


    char           integAxis  ;
    int      PMssFirstFileNum ;
    int      PMssLastFileNum  ;


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
    outputDir    = "";

      integAxis  = 'z';

       minMass    =  1e15;
       maxMass    = -1.0 ;
 radiusMultiplier =  1.0 ;
  radiusConverter = 1e-3 ;
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
     physicalSize = -1;

     particleMass = 0.0;
         redshift = -1;

//          setPartMass();

   maxIntegLength = 400.0;
      integStep   = 0.2;

        aexpn = 0;
       hubble = 0;
         om_o = 0;
         om_l = 0;

  PMssFirstFileNum = -1;
  PMssLastFileNum = -1;

}


class haloInfo {

  public:
//14 params

    void setX    ( float inpF ) {     x = inpF; }
    void setY    ( float inpF ) {     y = inpF; }
    void setZ    ( float inpF ) {     z = inpF; }
    void setC    ( float inpF ) {     C = inpF; }
    void setM    ( float inpF ) {     M = inpF; }
    void setN    ( float inpI ) {     N = inpI; }
    void setID   ( long  inpL ) {    id = inpL; }
    void setBA   ( float inpF ) { ba_rat= inpF; }
    void setCA   ( float inpF ) { ca_rat= inpF; }
    void setRm   ( float inpF ) { R_max = inpF; }
    void setPhi  ( float inpF ) {   phi = inpF; }
    void setTheta( float inpF ) {  theta= inpF; }

    void setDistinct( long inpL ) { distinct = inpL; }


    float getX    () const { return  x     ; }
    float getY    () const { return  y     ; }
    float getZ    () const { return  z     ; }
    float getC    () const { return  C     ; }
    float getM    () const { return  M     ; }
    float getN    () const { return  N     ; }
    float getRm   () const { return  R_max ; }
    long  getID   () const { return  id    ; }
    float getBA   () const { return  ba_rat; }
    float getCA   () const { return  ca_rat; }
    float getPhi  () const { return  phi   ; }
    float getTheta() const { return  theta ; }

    long  getDistinct() const { return  distinct; }


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

    float theta   ;
    float phi     ;

    // x horizonal axis, y verical, z along LOS
    // theta = angle oriented on xz,
    //         0 is oriented into z, pi/2 oriented into x
    //         range -pi/2 to pi/2
    // phi   = angle oriented on yz,
    //         0 is oriented into y, pi/2 oriented into z
    //         range 0 to pi


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


// Fortran function for reading PMss
extern "C" long long readpmss_( int  *jstep          , // Snapshot num
                                char *filestart      , // Part of the file name
                                int  *filestartlength, // Length of filestart
                                int  *firstPMss      , // First PMssFile Blah.0081.firstPMss.DAT
                                int  *lastPMss       );// Last ^


#endif // HE_CLASSES
