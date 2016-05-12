#include <cmath>
#include <vector>
#include <stdlib.h>


#include "halo_extraction_classes.h"
#include "triaxiality.h"



// Determines the orientation of the halo relative to our LOS
//   Writes orientation angles into the halo
int triaxiality(       haloInfo                  *halo ,   // The halo we are checking
                 const particlePosition     *particles ,   // The particle array
                 const long long           N_part_halo ,   // Number of particles we are using in the array
                 const long long          *partIndexes ) { // Indexes we are using in the array


  double I[3][3] = {  0 , 0 , 0 ,  // Inertia tensor
                      0 , 0 , 0 ,
                      0 , 0 , 0 };


  // Populate the inertia array with the particles
  for ( int i = 0; i < N_part_halo; ++i ){


  double m[3] = { (particles[  partIndexes[i]  ]).x_pos - (*halo).getX() ,  // Easier to represent particles as array when filling
                  (particles[  partIndexes[i]  ]).y_pos - (*halo).getY() ,  //  out the inertia tensor
                  (particles[  partIndexes[i]  ]).z_pos - (*halo).getZ() }; // Use indexes stored in partIndexes


  for ( int j = 0; j < 3          ; ++j ){
  for ( int k = 0; k < 3          ; ++k ){

    if ( j == k ){
      I[k][j] +=   ( m[0]*m[0]  + m[1]*m[1] + m[2]*m[2] - m[k]*m[j]  );
    }
    else{
      I[k][j] += - ( m[j]*m[k] );
    }
  }
  }
  }


  // The determinant for the eigenvalues evaluates to d L^3 + a L^2 + b L + c = 0
  //    with each constant being combinations of I

  double d = -1;

  double a = I[0][0]         + I[1][1]         + I[2][2];

  double b = I[0][1]*I[0][1] + I[0][2]*I[0][2] + I[1][2]*I[1][2] -
             I[0][0]*I[1][1] - I[0][0]*I[2][2] - I[1][1]*I[2][2] ;

  double c = I[0][0]*I[1][1]*I[2][2] + 2 * I[0][1]*I[0][2]*I[1][2] -
             I[0][0]*I[1][2]*I[1][2]   -   I[1][1]*I[0][2]*I[0][2] - I[2][2]*I[0][1]*I[0][1];


  // Zeros of this function must lie between and outside flat part of curve,
  //    take the derivative to find these locations with quadratic formula

  double zeroL, zeroR;

  {
    double A = 3 * d;
    double B = 2 * a;
    double C =     b;

    double zero1 = ( -B + sqrt( B*B - 4 * A * C ) ) / ( 2 * A );
    double zero2 = ( -B - sqrt( B*B - 4 * A * C ) ) / ( 2 * A );

    zeroL = std::min( zero1, zero2 );
    zeroR = std::max( zero1, zero2 );
  }

  double eigenValue[3];
  double eigenVector[3][3];

  // Use moving probes to find eigenvalues
  //   1st time through do left bound
  //   2nd time right bound
  //   3rd time center
  for ( int i = 0; i < 3; ++ i ){

    // Sets boundaries for probes based on region each run through
    double boundL, boundR;

    if ( 1 != setTriaxBounds( boundL, boundR, zeroL, zeroR, i, d, a, b, c ) ) {
      return -1;
    }


    double step;   // Step for the probes to move
    double diff;   // Differential to use for convergence
    int counter=0; // Counter for non-convergence


    // Step size will be dependent on distance between the two bounds
    step = ( boundR - boundL );

    // Move probes to find 0s
    do {

      step   =   step / 4.;                                // Step size needs to be reduced each runthrough
      boundL = moveProbe( boundL, boundR, step, d,a,b,c ); // Moves the left probe
      step   = - step;                                     // Flip sign on step for moving other probe
      boundR = moveProbe( boundR, boundL, step, d,a,b,c ); // Moves the right probe
      step   = - step;                                     // Return sign to normal
      diff   = boundR-boundL;                              // For checking boundary

      if ( diff < 0 ) counter = 100; // Went too far in moving probe

      ++counter;
    } while ((diff > 1e-8) && counter<100);

    if ( counter == 100 ){
      return -1;
    }

    eigenValue[i] = (boundR + boundL) / 2;

    double tempVector[3];

    calcEigenVector( tempVector, I, eigenValue[i] );

    for ( int j=0; j<3; ++j ){
      eigenVector[i][j] = tempVector[j];
    }

  } // 3 regions loop

  // Eigenvalues axes

  short  longAxisIndex(0);

  // Determine max/min eigenvales/axes
  if (   eigenValue[0] > eigenValue[1] ){
    if ( eigenValue[2] > eigenValue[0] ){ longAxisIndex = 1; } else
    if ( eigenValue[2] < eigenValue[1] ){ longAxisIndex = 2; } else{
                                          longAxisIndex = 1; }
  } else {
    if ( eigenValue[2] > eigenValue[1] ){ longAxisIndex = 0; } else
    if ( eigenValue[2] < eigenValue[0] ){ longAxisIndex = 2; } else{
                                          longAxisIndex = 0; }
  }


  // x horizonal axis, y verical, z along LOS
  // theta = angle oriented on xz,
  //         0 is oriented into z, pi/2 oriented into x
  //         range -pi/2 to pi/2
  // phi   = angle oriented on yz,
  //         0 is oriented into y, pi/2 oriented into z
  //         range 0 to pi

  double  x,  y,  z,  r;
  double  theta, phi;

  // Maxor axis direction
  x  = eigenVector[ longAxisIndex ][0];
  y  = eigenVector[ longAxisIndex ][1];
  z  = eigenVector[ longAxisIndex ][2];
  r  = std::sqrt( x*x + y*y + z*z );

  // Orientation angle
  theta = std::atan( x / z );
  phi   = std::acos( y / r );

  (*halo).setPhi(   phi   );
  (*halo).setTheta( theta );

  return 1;
}





// Calculates and returns eigenvector
void calcEigenVector(       double    eigenVector[3] ,  // Eigenvector to return
                      const double           i[3][3] ,  // Inertia array
                      const double           L       ){ // Eigenvalue to calc vectors for

  short swapAxis = -1;

  double I[3][3]; // Use a different I array since axes max need to be swapped

  // If any of the diagonal indexes are 0, make sure in the top left for easier calculations
  // 1 is x-y swap
  // 2 is x-z swap
  if ( i[1][1] == 0 ){       I[0][0] = i[1][1];       I[0][1] = i[1][0];       I[0][2] = i[1][2];
                             I[1][0] = i[0][1];       I[1][1] = i[0][0];       I[1][2] = i[0][2];
                             I[2][0] = i[2][1];       I[2][1] = i[2][0];       I[2][2] = i[2][2];
    swapAxis = 1;
  } else
  if ( i[2][2] == 0 ){       I[0][0] = i[2][2];       I[0][1] = i[2][1];       I[0][2] = i[2][0];
                             I[1][0] = i[1][2];       I[1][1] = i[1][1];       I[1][2] = i[1][0];
                             I[2][0] = i[0][2];       I[2][1] = i[0][1];       I[2][2] = i[0][0];
    swapAxis = 2;
  } else {                   I[0][0] = i[0][0];       I[0][1] = i[0][1];       I[0][2] = i[0][2];
                             I[1][0] = i[1][0];       I[1][1] = i[1][1];       I[1][2] = i[1][2];
                             I[2][0] = i[2][0];       I[2][1] = i[2][1];       I[2][2] = i[2][2];
  }


  double x, y, z;

  // More visually appealing array, with eigenvalues included

  double a1 = I[0][0]-L ;
  double b1 = I[0][1]   ;
  double b2 = I[1][1]-L ;
  double c1 = I[0][2]   ;
  double c2 = I[1][2]   ;
  double c3 = I[2][2]-L ;


  // Eigenvectors calculated based on zeros in array due to avoid dividing by 0,
  //   any zeros lead to special solutions,
  //   only physically real cases included
  if ( std::abs(a1) < 1e-8 ){

  if ( std::abs(b1) < 1e-8 &&
       std::abs(c3) < 1e-8 ){    z = 1;    y = - c2 / b2 * z;    x = - c2 / c1 * y;  } else
  if ( std::abs(c1) < 1e-8 ){    z = 1;    y = - c3 / c2 * z;    x = - c2 / b1 * z;  } else
  if ( std::abs(c2) < 1e-8 ){    z = 1;    y = - c1 / b1 * z;    x = - c3 / c1 * z;  } else
  if ( std::abs(c3) < 1e-8 ){    z = 1;    y = - c1 / b1 * z;    x = - c2 / c1 * y;  } else
                            {    z = 1;    y = - c1 / b1 * z;    x = - b2 / b1 * y
                                                                     - c2 / b1 * z; }

  } else
  if ( std::abs(b1) < 1e-8 ){    z = 1;    y = - c2 / b2 * z;    x = - c1 / a1 * z;  } else
  if ( std::abs(c1) < 1e-8 ){    z = 1;    y = - c3 / c2 * z;    x = - b1 / a1 * y;  } else
  if ( std::abs(c2) < 1e-8 ){    z = 1;    x = - c3 / c1 * z;    y = - b1 / b2 * x;  } else
  {

    // If no special case, use a long form solution
    z = 1;

    y = ( c3 / c1 - c2 / b1 ) / ( b2 / b1 - c2 / c1 ) * z;

    x = - b1 / a1 * y         -   c1 / a1 * z;

  }


  // Used to generate unit vectors
  double r = std::sqrt( x*x + y*y + z*z );

  // If we swapped values earlier when calculating the array,
  //   be sure to swap the eigenvector back
  if ( swapAxis == 1 ){

    eigenVector[0] = y / r;
    eigenVector[1] = x / r;
    eigenVector[2] = z / r;

  } else
  if ( swapAxis == 2 ){

    eigenVector[0] = z / r;
    eigenVector[1] = y / r;
    eigenVector[2] = x / r;

  }else{

    eigenVector[0] = x / r;
    eigenVector[1] = y / r;
    eigenVector[2] = z / r;

  }

  return;
}





// Find new position of mProbe,
//   moves new probe until sign matches stationary,
//   which indicates moved too far
double moveProbe( const double   mProbe ,  // Moving probe
                  const double   sProbe ,  // Stationary probe
                  const double     step ,  // Stepsize to move probe
                  const double        a ,  // Values of polynomial
                  const double        b ,
                  const double        c ,
                  const double        d ){


  double newProbe = mProbe;
  short   counter =      0;

  do{

    newProbe += step;

    if ( ++counter>5 ) return mProbe; // If for some reason went too far, abort
                                      //   only should be four steps max

  } while ( sign( cubicPoly( newProbe,a,b,c,d) ) !=
            sign( cubicPoly(   sProbe,a,b,c,d) ) );


  return newProbe-step;
}



// Returns sign of input
double sign( const double x ){
  if ( x > 0 ) return  1;
  if ( x < 0 ) return -1;
  return 0;
}



// Sets boundaries of zerofinding routine
// Returns 1 if set boundaries, -1 if couldn't
int setTriaxBounds(       double  &boundL ,   // Left  maximum boundary to return
                          double  &boundR ,   // Right maximum boundary to return
                    const double    zeroL ,   // Left  location of flat part of curve
                    const double    zeroR ,   // Right location of flat part of curve
                    const int      runNum ,   // 0-left zero, 1-right 0, else middle
                    const double        a ,   // Values of polynomial
                    const double        b ,
                    const double        c ,
                    const double        d ){


    double boundMod = zeroR-zeroL;

    // If zeros in same location, need to abort all
    if ( (1 - zeroL/zeroR) < 1e-6 ) return -1;

    if ( runNum == 0 ){ // Left zero

          int counter = 0;
               boundL = zeroL;
               boundR = zeroL;

          // Make sure we encompass the zero
          while ( true ){

            boundL = boundL - boundMod;

            if ( sign( cubicPoly( boundL,a,b,c,d) ) ==
                -sign( cubicPoly( boundR,a,b,c,d) ) ) break;

            // If not able to find a bound, exit
            if ( ++counter>10 ) return -1;
          }

    } else
    if ( runNum == 1 ){ // Right zero

          int counter = 0;
               boundL = zeroR;
               boundR = zeroR;

          // Make sure we encompass the zero
          while ( true ){

            boundR = boundR + boundMod;

            if ( sign( cubicPoly( boundL,a,b,c,d) ) ==
                -sign( cubicPoly( boundR,a,b,c,d) ) ) break;

            // If not able to find a bound, exit
            if ( ++counter>10 ) return -1;
          }


    } else{        // Center

               boundL = zeroL;
               boundR = zeroR;

          if ( sign( cubicPoly( boundL,a,b,c,d) ) ==
               sign( cubicPoly( boundR,a,b,c,d) ) ){
            return -1;
          }

    }

  return 1; // Successfully swapped
}


// Cubic polynomial
double cubicPoly( const double   x ,
                  const double   a ,
                  const double   b ,
                  const double   c ,
                  const double   d ){

  return a * x * x * x +
         b * x * x     +
         c * x         +
         d;

}

