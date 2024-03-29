#include <math.h>
#include <vector>

#include <CCfits/CCfits>

#include "halo_extraction_classes.h"
#include "link_halos.h"
#include "triaxiality.h"


//Sets the min/max values of the particle positions
void setMinMaxParticles( const particlePosition   *particle,  //Particles to find min/max of
                         inputInfo              &userInput ){ //Number of particles in the array

  //min and max values for each coordinate
  float minX( particle[0].x_pos );
  float maxX( particle[0].x_pos );
  float minY( particle[0].y_pos );
  float maxY( particle[0].y_pos );
  float minZ( particle[0].z_pos );
  float maxZ( particle[0].z_pos );

  //Find boundaries of particle positions
  for ( int i = 1; i < userInput.getNumParticles(); ++i ){
    minX = ( minX < particle[i].x_pos  ? minX : particle[i].x_pos);
    minY = ( minY < particle[i].y_pos  ? minY : particle[i].y_pos);
    minZ = ( minZ < particle[i].z_pos  ? minZ : particle[i].z_pos);
    maxX = ( maxX > particle[i].x_pos  ? maxX : particle[i].x_pos);
    maxY = ( maxY > particle[i].y_pos  ? maxY : particle[i].y_pos);
    maxZ = ( maxZ > particle[i].z_pos  ? maxZ : particle[i].z_pos);
  }

  //save the values
  userInput.setXmin( floor( minX ) );
  userInput.setYmin( floor( minY ) );
  userInput.setZmin( floor( minZ ) );
  userInput.setXmax(  ceil( maxX ) );
  userInput.setYmax(  ceil( maxY ) );
  userInput.setZmax(  ceil( maxZ ) );

}


// Locates halos that could be in our box with a full FOV
unsigned long findBoxHalos(       inputInfo  &userInput ,
                            const  haloInfo   *allHalos ,
                                   haloInfo   *boxHalos ,
                                      short    runNum   ){



  // Array to track how many halos are valid within the box
  short inBox[ userInput.getNumHalos() ];
  unsigned long sum(0);


  // Cycle through halos to determine if in our box
  for ( int i  = 0 ; i < userInput.getNumHalos(); ++i ){
      inBox[i] = 0;

    // If the halo is within the coordinates we have, and the FOV doesn't leave the box, save it
    if ( ( userInput.getXmin() < allHalos[i].getX() - userInput.getFOV() ) && ( allHalos[i].getX() + userInput.getFOV() < userInput.getXmax() ) &&
         ( userInput.getYmin() < allHalos[i].getY() - userInput.getFOV() ) && ( allHalos[i].getY() + userInput.getFOV() < userInput.getYmax() ) &&
         ( userInput.getZmin() < allHalos[i].getZ() - userInput.getFOV() ) && ( allHalos[i].getZ() + userInput.getFOV() < userInput.getZmax() )   ){

      inBox[i] = 1;
    }

    // Sum is the number of halos in the box
    sum += inBox[i];
  }


  // First time through, just finding number of valid halos
  if ( runNum == 0 ){
    return sum;
  }
  // Second time, do the saving
  else {

    // Go through and copy over the halos in the box to the new halo list
    unsigned long counter(0);
    for ( int i = 0 ; i < userInput.getNumHalos(); ++i ){
      if ( inBox[i] == 1 ){
        boxHalos[ counter ] = allHalos[ i ];
        ++counter;
      }
    }

    if ( sum != counter ){
      std::cout << " Error: Mismatch when selecting halos\n  Num in box:   " << sum << "\n  Num selected: " << counter << std::endl;

      logMessage( std::string( "Mismatch in counting halos in box: Num in box: "   ) +
                  std::to_string(  (long long) sum                                 ) +
                  std::string(     " Num selected: "                               ) +
                  std::to_string(  (long long) counter                             ) );

      logMessage( std::string("Aborting.") );

      exit(1);
    }

    userInput.setNumHalos( sum );


    return counter;
  }
}



//Make the link list of nearby particles
void makeLinkList( const inputInfo          userInfo ,   //Contains the global info needed
                   const particlePosition  *particle ,   //Array of the particles
                         long long         *myList   ,   //List pointing to next neighbor particle
                         long long         *myLabel  ){  //Label points to the last particle in list for index

  double cell = userInfo.getCell();

  int Nlx = userInfo.getNlx(); //Minimum box numbers (should be 0)
  int Nly = userInfo.getNly();
  int Nlz = userInfo.getNlz();
  int Nrx = userInfo.getNrx(); //Maximum box numbers
  int Nry = userInfo.getNry();
  int Nrz = userInfo.getNrz();


 //Fortran has -1, check if this matters, we are using unsigned long however
  for ( long i = 0 ; i < userInfo.getNumParticles() ; ++i ){
    myList[i] = -1;
  }

  for ( long i = 0 ; i < userInfo.getNtotCell() ; ++i ){
    myLabel[i] = -1;
  }


  //Loop over each particle, generating the link list
  for ( long i = 0; i < userInfo.getNumParticles() ; ++i ){

    //Find which box the particle is in, make sure we stay in the box
    int xIndex = std::min( std::max( Nlx, int( ( particle[i].x_pos - userInfo.getXmin() ) / cell ) ), Nrx );
    int yIndex = std::min( std::max( Nly, int( ( particle[i].y_pos - userInfo.getYmin() ) / cell ) ), Nry );
    int zIndex = std::min( std::max( Nlz, int( ( particle[i].z_pos - userInfo.getZmin() ) / cell ) ), Nrz );


    long lI = xIndex + yIndex * Nrx + zIndex * Nrx * Nry; //Label index


    myList[  i  ] = myLabel[ lI ]; //List previous particle in the same box
    myLabel[ lI ] = i;             //label stores last particle recorded in the box

  }

}



// Matches halos with particles, and calls write functions
void linkHaloParticles(              inputInfo   userInput ,  // Info from the user
                                      haloInfo      *halos ,  // The halo information
                        const particlePosition  *particles ,  // Position of the particles
                        const        long long  *labelList ,  // Label points to the first particle in cell
                        const        long long   *linkList ){ // Points to nearby neighbor particles


  int    Nlx = userInput.getNlx();  // Minimum box numbers (should be 0)
  int    Nly = userInput.getNly();
  int    Nlz = userInput.getNlz();
  int    Nrx = userInput.getNrx();  // Maximum box numbers
  int    Nry = userInput.getNry();
  int    Nrz = userInput.getNrz();

  float xMin = userInput.getXmin(); // Simplifying user values
  float yMin = userInput.getYmin();
  float zMin = userInput.getZmin();
  float cell = userInput.getCell();
  float FOV  = userInput.getFOV() ;

  double maxInteg        = userInput.getMaxIntegLength() / 2.0;                                     // Max distance to go on one side
  double    integStart   = userInput.getFOV() / 2.0;                                                // Starting value to integrate from, FOV/2
  int     N_integSteps   = ( log10( maxInteg ) - log10( integStart ) ) / userInput.getIntegStep() ; // Number of steps given our integration step size
  double    integStep    = ( log10( maxInteg ) - log10( integStart ) ) / N_integSteps;              // Refine the integration step size
  double   *integLengths = new double[ N_integSteps ];                                              // Array to fill with the different integration lengths


  logMessage( std::string( "Integration lengths:" ) );

  // Populate the array with the integration lengths
  for ( int i = 0 ; i < N_integSteps ; ++i ){
    integLengths[i] = pow( 10, (i+1) * integStep  +  log10(integStart) );

    logMessage( std::string(                  "   "            ) +
                std::to_string( (long double) integLengths[i]  ) );

  }



  // If we have the whole z axis, use periodic boundary conditions
  short periodic = 0;
  if ( userInput.getZmax() - zMin == userInput.getCatBoxSize() ){
        periodic = 1;
    std::cout<<"  Using periodic boundaries..." << std::endl;

    logMessage( std::string( "Using periodic boundaries" ) );

  }


  // Cycle through each halo, trying to locate the needed particles
  for ( int i = 0; i < userInput.getNumHalos(); ++i ){

    // Simplify coordinates
    float x   = halos[i].getX();
    float y   = halos[i].getY();
    float z   = halos[i].getZ();
    float r_v = halos[i].getRm() * userInput.getRadiusConvert();

    // Left and rightmost coordinates in the FOV
    float xLeft  = x - FOV / 2.0;
    float yLeft  = y - FOV / 2.0;
    float zLeft  = z - FOV / 2.0;
    float xRight = x + FOV / 2.0;
    float yRight = y + FOV / 2.0;
    float zRight = z + FOV / 2.0;


    // Box to the left and right of the box the center of the halo is located in
    int xLeftN  = std::min( std::max( Nlx, int( ( x - xMin            ) / cell - 1 ) ), Nrx );
    int yLeftN  = std::min( std::max( Nly, int( ( y - yMin            ) / cell - 1 ) ), Nry );
    int zLeftN  = std::min( std::max( Nlz, int( ( z - zMin - maxInteg ) / cell - 1 ) ), Nrz ); // z's need to include the integration lengths
    int xRightN = std::min( std::max( Nlx, int( ( x - xMin            ) / cell + 1 ) ), Nrx );
    int yRightN = std::min( std::max( Nly, int( ( y - yMin            ) / cell + 1 ) ), Nry );
    int zRightN = std::min( std::max( Nlz, int( ( z - zMin + maxInteg ) / cell + 1 ) ), Nrz );


    if ( periodic==1 ){
                        zLeftN =           int( ( z - zMin - maxInteg ) / cell - 1 )         ; // If loop out of bounds, it's alright
                       zRightN =           int( ( z - zMin + maxInteg ) / cell + 1 )         ;

      if ( zLeftN < 0 ) {
          zLeftN =   zLeftN + Nrz + 1;  // If left side negative, need to add n nodes to each side
         zRightN =  zRightN + Nrz + 1;
      }
    }


    // Number of particles in each set
    long long    N_sphere(0);
    long long    N_box   (0);
    long long   *N_integ  =  new long long [ N_integSteps ]; // N_integ will contain the number in each integration set
    long long maxN_integ  = 0;                               // maxN_integ will contain the max index


    // Maximum possible integration length for our halo
    double haloMaxInteg      = maxInteg;
    int    haloMaxIntegIndex = N_integSteps;



    // Find if there is an integ length that goes outside the box, keep looping until it's in
    for ( int j = N_integSteps-1; j >= 0 ; --j ){
      N_integ[j] = 0;

      if ( ( periodic !=1 )                              && (
           ( z + integLengths[j] ) > userInput.getZmax() || // If an integration length is outside the box, and not periodic
           ( z - integLengths[j] ) < userInput.getZmin() )){
        if ( j > 0 ){
          haloMaxInteg      = integLengths[j-1];
          haloMaxIntegIndex = j-1;
        } else {
          haloMaxInteg      = integStart;
          haloMaxIntegIndex = -1;
        }
      }
    }



    // Loop over cells possibly containing particles of our halo, and count
    //  the number that belong in each particle set
    for ( int xIndex = xLeftN; xIndex <= xRightN; ++xIndex ){
    for ( int yIndex = yLeftN; yIndex <= yRightN; ++yIndex ){
    for ( int zI     = zLeftN; zI     <= zRightN; ++zI     ){

      int zIndex = zI % (Nrz+1);                         // Handles boundaries

      long cellIndex = xIndex +                          //Index for the cell for labellist
                       yIndex * Nrx +
                       zIndex * Nrx * Nry;

      long long particleIndex = labelList[ cellIndex ];  //The index for the first particle in the link list

      //Loop over link list
      while ( particleIndex != -1 ){

        float partX = particles[ particleIndex ].x_pos;
        float partY = particles[ particleIndex ].y_pos;
        float partZ = particles[ particleIndex ].z_pos;


        // Check if particle is in a set

        if ( (  xLeft              < partX) && (partX < ( xRight           )) &&
             (  yLeft              < partY) && (partY < ( yRight           )) &&
             (( z - haloMaxInteg ) < partZ) && (partZ < ( z + haloMaxInteg )) ){

          // Check if in sphere, Rm > sqrt( ( x_h-x_p )^2 + ...
          if ( ( r_v * r_v )  >= ( ( x - partX ) * ( x - partX ) +
                                   ( y - partY ) * ( y - partY ) +
                                   ( z - partZ ) * ( z - partZ ) ) ){
            ++N_sphere;
          }

          // Particle in the FOV frame
          else if ( ( zLeft < partZ ) && ( partZ < zRight ) ) {
            ++N_box;
          }

          // Find the integration length it is in
          else {
            int integIndex = std::min( std::max(    int( ( log10( abs(z - partZ) ) - log10(integStart) ) / integStep )     , 0 ) , N_integSteps - 1 );
            N_integ[ integIndex ] += 1;

            if ( N_integ[integIndex] > maxN_integ ){
              maxN_integ = N_integ[ integIndex ];       // Find the maximum amount of particles in any integ set
            }
          }
        }

        particleIndex = linkList[ particleIndex ];
      } //while over link list
    }
    }   //Box loop
    }

    // Only continue if we found particles for our halo
    if ( ( N_sphere > 0 ) && ( N_box > 0 ) ){


      logMessage( std::string(    "Halo: "         ) +
                  std::to_string( (long long) halos[i].getID() ) +
                  std::string(    ", N_sphere= "   ) +
                  std::to_string(    N_sphere      ) +
                  std::string(    ", N_box= "      ) +
                  std::to_string(    N_box         ) );


      // Need to allocate for storing indexes
      long long *sphereIndexes = new long long [ N_sphere ];
      long long *   boxIndexes = new long long [ N_box    ];
      long long * integIndexes = new long long [ N_integSteps * maxN_integ ];

      for ( int kk = 0; kk < N_sphere                  ; ++kk )
        sphereIndexes[kk] = 0;
      for ( int kk = 0; kk < N_box                     ; ++kk )
           boxIndexes[kk] = 0;
      for ( int kk = 0; kk < N_integSteps * maxN_integ ; ++kk )
         integIndexes[kk] = 0;


      //Counter for the index of each
      long long      sCounter(0);
      long long      bCounter(0);
      long long *integCounter = new long long [ N_integSteps ]; // N_integ just contains the number in each integration set


      for ( int foo = 0; foo < N_integSteps; ++foo ){
        integCounter[foo] = 0;
      }

      //Loop over cells again to find indexes of particles
      for ( int xIndex = xLeftN; xIndex <= xRightN; ++xIndex ){
      for ( int yIndex = yLeftN; yIndex <= yRightN; ++yIndex ){
      for ( int zI     = zLeftN; zI     <= zRightN; ++zI     ){

        int zIndex = zI % (Nrz+1);                         // Handles periodic boundaries

        long cellIndex = xIndex +                          //Index for the cell for labellist
                         yIndex * Nrx +
                         zIndex * Nrx * Nry;

        long long particleIndex = labelList[ cellIndex ];  //The index for the first particle in the link list


        //Loop over link list
        while ( particleIndex != -1 ){


          float partX = particles[ particleIndex ].x_pos;
          float partY = particles[ particleIndex ].y_pos;
          float partZ = particles[ particleIndex ].z_pos;

          if ( (( xLeft            ) < partX) && (partX < ( xRight           )) &&
               (( yLeft            ) < partY) && (partY < ( yRight           )) &&
               (( z - haloMaxInteg ) < partZ) && (partZ < ( z + haloMaxInteg )) ){


            // Check if in sphere, Rm > sqrt( ( x_h-x_p )^2 + ...
            if ( ( r_v * r_v ) >= ( ( x - partX ) * ( x - partX ) +
                                    ( y - partY ) * ( y - partY ) +
                                    ( z - partZ ) * ( z - partZ ) ) ){
              sphereIndexes[ sCounter ] = particleIndex;
              ++sCounter;
            }

            // Particle in the FOV frame
            else if ( ( zLeft < partZ ) && ( partZ < zRight ) ) {

              boxIndexes[ bCounter ] = particleIndex;
              ++bCounter;
            }

            // Find the integration length it is in
            else {

              int integIndex = std::min( std::max(    int( ( log10( abs(z - partZ) ) - log10(integStart) ) / integStep )     , 0 ) , N_integSteps - 1 );

              integIndexes[ integCounter[ integIndex ] + integIndex * maxN_integ ] = particleIndex;

              integCounter[ integIndex ] += 1;
            }
          }
          particleIndex = linkList[ particleIndex ];
        } //while
      }
      }   //Box loops
      }


      //Write the fits file for this halo
      writeImage( userInput       ,
                  &halos[i]       ,
                  particles       ,
                   N_sphere       ,
                   N_box          ,
                   N_integ        ,
                maxN_integ        ,
                   N_integSteps   ,
                     integLengths ,
                    sphereIndexes ,
                       boxIndexes ,
                     integIndexes );

    delete[] sphereIndexes, boxIndexes, integIndexes, integCounter;
    }     //If we found particles
    delete [] N_integ;

  }       //Halo loop

  delete [] integLengths;

}



// Attempts to write the many images
void writeImage(              inputInfo    userInput  , // All the user info
                               haloInfo        *halo  , // The halo we are considering
                 const particlePosition  *particles   , // The full array of particles
                 const long long          N_sphere    ,
                 const long long          N_box       , // Number of particles in sphere set, box set, and integration length sets
                 const long long         *N_integ     ,
                 const long long       maxN_integ     , // Maximum number on particles in any integIndexes bin
                 const long long          N_integBins , // Number of bins (steps) of integration
                 const    double      *integLengths   , // Contains the integration lengths
                 const long long     *sphereIndexes   , // Particle indexes in each set
                 const long long     *   boxIndexes   ,
                 const long long     * integIndexes   ){



  // Number of elements in each direction, and total number of pixels
  int   N_pixels[2] = { userInput.getNPixlesH(), userInput.getNPixelsV() };
  long lN_pixels[2] = { N_pixels[0], N_pixels[1] };

  int   N_elements  = std::accumulate(&N_pixels[0],&N_pixels[2],1,std::multiplies<int>());


  // Our surface density array
  std::valarray<double> SD( N_elements );
  for ( int i = 0; i < N_elements; ++i ){
    SD[i] = 0;
  }

  // Determines triaxiality and orientation angles along LOS
  //   halo saves the angle info
  // If couldn't find info, return obvious flag
  if ( triaxiality( halo, particles, N_sphere, sphereIndexes ) == -1 ) {
    (*halo).setPhi(   -100 );
    (*halo).setTheta( -100 );

    logMessage( std::string(    "Unable to determine primary axes for halo: ") +
                std::to_string( (long long)                  (*halo).getID() ) );
  }

  // Generate file names, require temp for the sprintf
  char temp[100];


  // Write the halo file name
  sprintf( temp,              "%sHalo%s_%09li_%04.1f.FITS", (userInput.getParticleDir()).c_str(), (userInput.getCatType()).c_str(), (*halo).getID(), userInput.getFOV() );
  const std::string sphereFileName = temp;



  // Attempt to write the file, if it fails abort
  if (  writeFits( sphereFileName ,
                      N_pixels    ,
                      N_elements  ,
                              &SD ,
                      N_sphere    ,
                    sphereIndexes ,
                        particles ,
                            *halo ,
                        userInput ,
                                0 ) == -1 ){

    std::cout<<"Failed to write file: "<<sphereFileName<<std::endl;
    exit(1);
  }
    std::cout<<"    Wrote file: "      <<sphereFileName<<std::endl;


  // Write the box file name
  sprintf( temp,        "%sBox%s_%09li_%04.1f.FITS", (userInput.getParticleDir()).c_str(), (userInput.getCatType()).c_str(), (*halo).getID(), userInput.getFOV() );
  const std::string    boxFileName = temp;

  // Attempt to write the file, if it fails abort
  if (  writeFits(    boxFileName ,
                      N_pixels    ,
                      N_elements  ,
                              &SD ,
                      N_box       ,
                       boxIndexes ,
                        particles ,
                            *halo ,
                        userInput ,
              userInput.getFOV()  ) == -1 ){

    std::cout<<"Failed to write file: "<<   boxFileName<<std::endl;
    exit(1);
  }
    std::cout<<"    Wrote file: "      <<   boxFileName<<std::endl;

  // Loop over integration lengths, if there are particles in the bin, write that index file
  if ( maxN_integ    > 0 ){
  for ( int i = 0; i < N_integBins; ++i ){
  if (    N_integ[i] > 0 ){

    // Integration file names
    sprintf( temp, "%sBox%s_%09li_%04.1f_%06.1lf.FITS",
        (userInput.getParticleDir()).c_str(), (userInput.getCatType()).c_str(), (*halo).getID(), userInput.getFOV(), 2.0* integLengths[i] );
    const std::string iFileName = temp;


    // Need a 1d array for writefiles, so only take those indexes
    long long iIndexes[ 2 ];
    iIndexes[0] = 0;

    float x_start = (*halo).getX() - userInput.getFOV() / 2.0;
    float y_start = (*halo).getY() - userInput.getFOV() / 2.0;

    // Scaling of pixels on image
    float x_scale = N_pixels[0] / userInput.getFOV() ;
    float y_scale = N_pixels[1] / userInput.getFOV() ;


    // Place our particles in the SD array before entering function
    // This is due to segfaulting the memory for large FOVs
    for ( int ii = 0; ii < N_integ[i]; ++ii ){


      // Index of particle in the set
      long long index = integIndexes[ii];

      // Location in the image, physical coordinates
      float x_pos = particles[index].x_pos - x_start;
      float y_pos = particles[index].y_pos - y_start;

      // Which pixel it's in
      int  hIndex = std::min( std::max(  int( x_pos * x_scale ) , 0 ), N_pixels[0] );
      int  vIndex = std::min( std::max(  int( y_pos * y_scale ) , 0 ), N_pixels[1] );

      // Place it
      SD[ vIndex * N_pixels[0] + hIndex ] += userInput.getParticleMass();
    }


    // Attempt to write the file, if it fails abort
    if (  writeFits(      iFileName ,
                        N_pixels    ,
                        N_elements  ,
                                &SD ,
                                0   ,
                           iIndexes ,
                          particles ,
                              *halo ,
                          userInput ,
              2.0 * integLengths[i] ) == -1 ){

    std::cout<<"Failed to write file: "<<     iFileName<<std::endl;
    exit(1);
  }
    std::cout<<"    Wrote file: "      <<     iFileName<<std::endl;

  }
  }
  }
    std::cout << std::endl;


}


// Write the fits image to file
int  writeFits( const std::string           fileName    ,  // File name to write
                const int                  *N_pixels    ,  // N_pixels in each direction
                const int                   N_pixelsTot ,  // Total number of pixels
                std::valarray<double>             *SD   ,  // SD array to use, passed because we will keep adding to the array
                const long long             N_indexes   ,  // Number of indexes in set to add
                const long long              *indexes   ,  // Indexes of the set
                const particlePosition     *particles   ,  // All the particles
                const         haloInfo           halo   ,  // Central halo
                             inputInfo      userInput   ,  // All the user information
                const float               integLength   ){ // Integration length, if applicable


  // Needs to be long format for functions
  long lN_pixels[2] = { N_pixels[0], N_pixels[1] };


  // Furthest left locations on image in real space
  float x_start = halo.getX() - userInput.getFOV() / 2.0;
  float y_start = halo.getY() - userInput.getFOV() / 2.0;

  // Scaling of pixels on image
  float x_scale = N_pixels[0] / userInput.getFOV() ;
  float y_scale = N_pixels[1] / userInput.getFOV() ;


  // Place our particles in the SD array
  for ( int i = 0; i < N_indexes; ++i ){

    // Index of particle in the set
    long long index = indexes[i];

    // Location in the image, physical coordinates
    float x_pos = particles[index].x_pos - x_start;
    float y_pos = particles[index].y_pos - y_start;

    // Which pixel it's in
    int  hIndex = std::min( std::max(  int( x_pos * x_scale ) , 0 ), N_pixels[0] );
    int  vIndex = std::min( std::max(  int( y_pos * y_scale ) , 0 ), N_pixels[1] );

    // Place it
    (*SD)[ vIndex * N_pixels[0] + hIndex ] += userInput.getParticleMass();
  }


  // Get total mass in image
  float totM = 0;
  for ( int i = 0; i < N_pixelsTot; ++i ){
    totM += (*SD)[i];
  }


  using namespace CCfits;
  std::auto_ptr<FITS> pFits(0);

  try
  {
    const std::string afileName( "!"+fileName );
    pFits.reset( new FITS( afileName, DOUBLE_IMG, 2, lN_pixels) );
  }
  catch (FITS::CantCreate){  return -1; }

  int axisNum =0;

  if ( userInput.getIntegAxis() == 'x' ) axisNum = 1;
  if ( userInput.getIntegAxis() == 'y' ) axisNum = 2;
  if ( userInput.getIntegAxis() == 'z' ) axisNum = 3;

  // Add a header
  ( *pFits ).pHDU().addKey("Catalog"     , userInput.getCatType()      , "Simulation catalog");
  ( *pFits ).pHDU().addKey("FOV"         , userInput.getFOV()          , "Field of view, (h^{-1} Mpc)");
  ( *pFits ).pHDU().addKey("N_pixels_v"  , userInput.getNPixelsV()     , "Number of pixels on y axis");
  ( *pFits ).pHDU().addKey("N_pixels_h"  , userInput.getNPixlesH()     , "Number of pixels on x axis");
  ( *pFits ).pHDU().addKey("Integ axis"  , axisNum                     , "Axis parallel to LOS");
  ( *pFits ).pHDU().addKey("TotMass"     , totM                        , "Total mass in the image");
  ( *pFits ).pHDU().addKey("PHYSICALSIZE", userInput.getPhysicalSize() , "Angular size designated by the user");
  ( *pFits ).pHDU().addKey("PIXELUNIT"   , 1./userInput.getPixelUnit() , "Mass unit of pixels");
  ( *pFits ).pHDU().addKey("REDSHIFT"    , userInput.getRedshift()     , "Redshift designated by the user");


  ( *pFits ).pHDU().addKey("X"       , halo.getX()      , "X coordinate of the halo, center of image");
  ( *pFits ).pHDU().addKey("Y"       , halo.getY()      , "Y coordinate of the halo, center of image");
  ( *pFits ).pHDU().addKey("Z"       , halo.getZ()      , "Z coordinate of the halo, center of image");
  ( *pFits ).pHDU().addKey("ID"      , halo.getID()     ,          "ID number of the central halo");
  ( *pFits ).pHDU().addKey("Mass"    , halo.getM()      ,               "Mass of the central halo (M_{\\odot})");
  ( *pFits ).pHDU().addKey("Rvir"    , halo.getRm()     ,      "Virial radius of the central halo (h^{-1}kpc)");
  ( *pFits ).pHDU().addKey("C"       , halo.getC()      ,      "Virial radius of the central halo (h^{-1}kpc)");
  ( *pFits ).pHDU().addKey("b/a"     , halo.getBA()     ,   "Ratio of b to a axis of central halo");
  ( *pFits ).pHDU().addKey("c/a"     , halo.getCA()     ,   "Ratio of c to a axis of central halo");
  ( *pFits ).pHDU().addKey("phi"     , halo.getPhi()    ,   "Angle on yz plane, 0 along +y, pi/2 along +z");
  ( *pFits ).pHDU().addKey("theta"   , halo.getTheta()  ,   "Angle on xz plane, 0 along +z, pi/2 along +x");


  // When doing this with an integration length, add to header
  if ( integLength > 0 )
  ( *pFits ).pHDU().addKey("Integ", integLength , "Integration length along LOS");

  // Write the image
  ( *pFits ).pHDU().write( 1, N_pixelsTot, *SD);


  logMessage( std::string( "    Wrote file: " ) +
                                fileName        );


  userInput.wroteFile();
  return 0;
}



