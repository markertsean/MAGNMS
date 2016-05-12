
#include <halo_extraction_classes.h>

#ifndef TRIAXIALITY
#define TRIAXIALITY


int triaxiality(       haloInfo                  *halo ,   // The halo we are checking
                 const particlePosition     *particles ,   // The particle array
                 const long long           N_part_halo ,   // Number of particles we are using in the array
                 const long long          *partIndexes );  // Indexes we are using in the array


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



#endif // TRIAXIALITY
