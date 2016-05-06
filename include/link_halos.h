
#include "halo_extraction_classes.h"

#ifndef LINK_HALOS
#define LINK_HALOS



//Sets the min/max values of the particle positions
void setMinMaxParticles( const particlePosition  *particle ,  //Particles to find min/max of
                         inputInfo              &userInput ); //Number of particles in the array

//Locates halos that could be in our box with a full FOV
unsigned long findBoxHalos(      inputInfo  &userInfo  ,
                            const haloInfo  *allHalos  ,
                                  haloInfo  *boxHalos  ,
                                     short   runNum    );

//Make the link list of nearby particles
void makeLinkList( const inputInfo          userInfo ,   //Contains the global info needed
                   const particlePosition  *particle ,   //Array of the particles
                         long long         *myList   ,   //List pointing to next neighbor particle
                         long long         *myLabel  );  //Label points to the last particle in list for index


void linkHaloParticles(              inputInfo   userInput ,  // Info from the user
                        const         haloInfo      *halos ,  // The halo information
                        const particlePosition  *particles ,  // Position of the particles
                        const        long long  *labelList ,  // Label points to the first particle in cell
                        const        long long   *linkList ); // Points to nearby neighbor particles


void writeImage(              inputInfo    userInput  , // All the user info
                 const         haloInfo         halo  , // The halo we are considering
                 const particlePosition  *particles   , // The full array of particles
                 const long long          N_sphere    ,
                 const long long          N_box       , // Number of particles in sphere set, box set, and integration length sets
                 const long long         *N_integ     ,
                 const long long       maxN_integ     , // Maximum number on particles in any integIndexes bin
                 const long long          N_integBins , // Number of bins (steps) of integration
                 const    double      *integLengths   , // Contains the integration lengths
                 const long long     *sphereIndexes   , // Particle indexes in each set
                 const long long     *   boxIndexes   ,
                 const long long     * integIndexes   );

// Write the fits image to file
int  writeFits( const std::string           fileName            ,  // File name to write
                const int                  *N_pixels            ,  // N_pixels in each direction
                const int                   N_pixelsTot         ,  // Total number of pixels
                std::valarray<double>             *SD           ,  // SD array to use, passed because we will keep adding to the array
                const long long             N_indexes           ,  // Number of indexes in set to add
                const long long              *indexes           ,  // Indexes of the set
                const particlePosition     *particles           ,  // All the particles
                const         haloInfo           halo           ,  // Central halo
                             inputInfo      userInput           ,  // All the user information
                const float               integLength = -1.0    ); // Integration length, if applicable



int triaxiality();

// Returns the sign of a variable as 1, -1, or 0
double sign( const double x );

// Sets boundaries of zerofinding routine
// Returns 1 if set boundaries, -1 if couldn't
int setTriaxBounds( double  &boundL ,   // Left  maximum boundary to return
                    double  &boundR ,   // Right maximum boundary to return
              const double    zeroL ,   // Left  location of flat part of curve
              const double    zeroR ,   // Right location of flat part of curve
              const int      runNum ,   // 0-left zero, 1-right 0, else middle
              const double        a ,   // Values of polynomial
              const double        b ,
              const double        c ,
              const double        d );


// Cubic polynomial
double cubicPoly( const double   x ,
                  const double   a ,
                  const double   b ,
                  const double   c ,
                  const double   d );

// Find new position of mProbe,
//   moves new probe until sign matches stationary,
//   which indicates moved too far
double moveProbe( const double   mProbe ,  // Moving probe
                  const double   sProbe ,  // Stationary probe
                  const double     step ,  // Stepsize to move probe
                  const double        a ,  // Values of polynomial
                  const double        b ,
                  const double        c ,
                  const double        d );

void calcEigenVector(       double    eigenVector[3] ,  // Eigenvector to return
                      const double           i[3][3] ,  // Inertia array
                      const double           L       ); // Eigenvalue to calc vectors for



#endif // LINK_HALOS
