#include <math.h>

#include <CCfits/CCfits>

#include "halo_extraction_classes.h"
#include "link_halos.h"


//Sets the min/max values of the particle positions
void setMinMaxParticles( particlePosition particle[],  //Particles to find min/max of
                         inputInfo       &userInput ){ //Number of particles in the array

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




//Locates halos that could be in our box with a full FOV
unsigned long findBoxHalos( inputInfo &userInput  ,
                             haloInfo  allHalos[] ,
                             haloInfo  boxHalos[] ,
                               short   runNum     ){



  //Array to track how many halos are valid within the box
  short inBox[ userInput.getNumHalos() ];
  unsigned long sum(0);


  //Cycle through halos to determine if in our box
  for ( int i  = 0 ; i < userInput.getNumHalos(); ++i ){
      inBox[i] = 0;

    //If the halo is within the coordinates we have, and the FOV doesn't leave the box, save it
    if ( ( userInput.getXmin() < allHalos[i].getX() - userInput.getFOV() ) && ( allHalos[i].getX() + userInput.getFOV() < userInput.getXmax() ) &&
         ( userInput.getYmin() < allHalos[i].getY() - userInput.getFOV() ) && ( allHalos[i].getY() + userInput.getFOV() < userInput.getYmax() ) &&
         ( userInput.getZmin() < allHalos[i].getZ() - userInput.getFOV() ) && ( allHalos[i].getZ() + userInput.getFOV() < userInput.getZmax() )   ){

      inBox[i] = 1;
    }

    //Sum is the number of halos in the box
    sum += inBox[i];
  }


  //First time through, just finding number of valid halos
  if ( runNum == 0 ){
    return sum;
  }
  //Second time, do the saving
  else {

    //Go through and copy over the halos in the box to the new halo list
    unsigned long counter(0);
    for ( int i = 0 ; i < userInput.getNumHalos(); ++i ){
      if ( inBox[i] == 1 ){
        boxHalos[ counter ] = allHalos[ i ];
        ++counter;
      }
    }

    if ( sum != counter ){
      std::cout << " Error: Mismatch when selecting halos\n  Num in box:   " << sum << "\n  Num selected: " << counter << std::endl;
      exit(1);
    }

    userInput.setNumHalos( sum );


    return counter;
  }
}



//Make the link list of nearby particles
void makeLinkList( inputInfo        userInfo   ,   //Contains the global info needed
                   particlePosition particle[] ,   //Array of the particles
                   long long        myList  [] ,   //List pointing to next neighbor particle
                   long long        myLabel [] ){  //Label points to the last particle in list for index

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

  for ( int i = Nlx; i < Nrx ; ++i ){  //x
  for ( int j = Nly; j < Nry ; ++j ){  //y
  for ( int k = Nlz; k < Nrz ; ++k ){  //z
    myLabel[ i + j * Nrx + k * Nrx * Nry ] = -1;
  }
  }
  }

  //Loop over each particle, generating the link list
  for ( long i = 0; i < 100;++i){//userInfo.getNumParticles() ; ++i ){

    //Find which box the particle is in, make sure we stay in the box
    int xIndex = std::min( std::max( Nlx, int( ( particle[i].x_pos - userInfo.getXmin() ) / cell ) ), Nrx );
    int yIndex = std::min( std::max( Nly, int( ( particle[i].y_pos - userInfo.getYmin() ) / cell ) ), Nry );
    int zIndex = std::min( std::max( Nlz, int( ( particle[i].z_pos - userInfo.getZmin() ) / cell ) ), Nrz );


    long lI = xIndex + yIndex * Nrx + zIndex * Nrx * Nry; //Label index


    myList[  i  ] = myLabel[ lI ]; //List previous particle in the same box
    myLabel[ lI ] = i;             //label stores last particle recorded in the box
//printf("%5li %5li %5li %5li\n",i,lI,myList[i],myLabel[lI]);

  }

}



//Matches halos with particles, and calls write functions
void linkHaloParticles( inputInfo userInput   ,
                         haloInfo     halos[] ,
                 particlePosition particles[] ,
                       long long  labelList[] ,
                       long long   linkList[] ){




  int Nlx = userInput.getNlx(); //Minimum box numbers (should be 0)
  int Nly = userInput.getNly();
  int Nlz = userInput.getNlz();
  int Nrx = userInput.getNrx(); //Maximum box numbers
  int Nry = userInput.getNry();
  int Nrz = userInput.getNrz();

  float xMin = userInput.getXmin();
  float yMin = userInput.getYmin();
  float zMin = userInput.getZmin();
  float cell = userInput.getCell();
  float FOV  = userInput.getFOV() ;




  for ( int i = 0; i < userInput.getNumHalos(); ++i ){




    float x = halos[i].getX();
    float y = halos[i].getY();
    float z = halos[i].getZ();


    //Box to the left and right of the box the center of the halo is located in
    int xLeft  = std::min( std::max( Nlx, int( ( halos[i].getX() - xMin ) / cell - 1 ) ), Nrx );
    int yLeft  = std::min( std::max( Nly, int( ( halos[i].getY() - yMin ) / cell - 1 ) ), Nry );
    int zLeft  = std::min( std::max( Nlz, int( ( halos[i].getZ() - zMin ) / cell - 1 ) ), Nrz );
    int xRight = std::min( std::max( Nlx, int( ( halos[i].getX() - xMin ) / cell + 1 ) ), Nrx );
    int yRight = std::min( std::max( Nly, int( ( halos[i].getY() - yMin ) / cell + 1 ) ), Nry );
    int zRight = std::min( std::max( Nlz, int( ( halos[i].getZ() - zMin ) / cell + 1 ) ), Nrz );


    //Number of particles in each set
    unsigned long N_sphere(0);
    unsigned long N_box   (0);
    unsigned long N_i1    (0);
    unsigned long N_i2    (0);
    unsigned long N_i3    (0);


    //Loop over cells possibly containing particles of our halo, and count
    // the number that belong in each particle set
    for ( int xIndex = xLeft; xIndex <= xRight; ++xIndex ){
    for ( int yIndex = yLeft; yIndex <= yRight; ++yIndex ){
    for ( int zIndex = zLeft; zIndex <= zRight; ++zIndex ){



      long cellIndex = xIndex +                          //Index for the cell for labellist
                       yIndex * Nrx +
                       zIndex * Nrx * Nry;

      long long particleIndex = labelList[ cellIndex ];  //The index for the first particle in the link list



      //Loop over link list
      while ( particleIndex != -1 ){


        float partX = particles[ particleIndex ].x_pos;
        float partY = particles[ particleIndex ].y_pos;
        float partZ = particles[ particleIndex ].z_pos;


        //Check if in sphere, Rm > sqrt( ( x_h-x_p )^2 + ...
        if ( ( halos[i].getRm() * halos[i].getRm() ) > ( ( halos[i].getX() - partX ) * ( halos[i].getX() - partX ) +
                                                         ( halos[i].getY() - partY ) * ( halos[i].getY() - partY ) +
                                                         ( halos[i].getZ() - partZ ) * ( halos[i].getZ() - partZ ) ) ){
          ++N_sphere;
        }

        //Particle in the box frame
        else if ( ( (halos[i].getX()-FOV) < partX ) && ( partX < (halos[i].getX()+FOV) ) &&
                  ( (halos[i].getY()-FOV) < partY ) && ( partY < (halos[i].getY()+FOV) ) &&
                  ( (halos[i].getZ()-FOV) < partZ ) && ( partZ < (halos[i].getZ()+FOV) ) ){
          ++N_box;
        }
//Integration length, needs modification of z index



        particleIndex = linkList[ particleIndex ];
      } //while
    }
    }   //Box loop
    }


//N_box = 1;
    //Only continue if we found particles for our halo
    if ( ( N_sphere > 0 ) && ( N_box > 0 ) ){
//N_box = 0;


      //Need to allocate for indexes
      long long *sphereIndexes = new long long [ N_sphere ];
      long long *   boxIndexes = new long long [ N_box    ];
      long long *    i1Indexes;
      long long *    i2Indexes;
      long long *    i3Indexes;

      //If we were able to complete any integration, use it
      if ( N_i1 > 0 ) i1Indexes = new long long [ N_i1 ];
      if ( N_i2 > 0 ) i2Indexes = new long long [ N_i2 ];
      if ( N_i3 > 0 ) i3Indexes = new long long [ N_i3 ];

      //Counter for the index of each
      long long  sCounter(0);
      long long  bCounter(0);
      long long i1Counter(0);
      long long i2Counter(0);
      long long i3Counter(0);


      //Loop over cells again to find indexes of particles
      for ( int xIndex = xLeft; xIndex <= xRight; ++xIndex ){
      for ( int yIndex = yLeft; yIndex <= yRight; ++yIndex ){
      for ( int zIndex = zLeft; zIndex <= zRight; ++zIndex ){



        long cellIndex = xIndex +                          //Index for the cell for labellist
                         yIndex * Nrx +
                         zIndex * Nrx * Nry;

        long long particleIndex = labelList[ cellIndex ];  //The index for the first particle in the link list


        //Loop over link list
        while ( particleIndex != -1 ){


          float partX = particles[ particleIndex ].x_pos;
          float partY = particles[ particleIndex ].y_pos;
          float partZ = particles[ particleIndex ].z_pos;


          //Check if in sphere, Rm > sqrt( ( x_h-x_p )^2 + ...
          if ( ( halos[i].getRm() * halos[i].getRm() ) > ( ( halos[i].getX() - partX ) * ( halos[i].getX() - partX ) +
                                                           ( halos[i].getY() - partY ) * ( halos[i].getY() - partY ) +
                                                           ( halos[i].getZ() - partZ ) * ( halos[i].getZ() - partZ ) ) ){
            sphereIndexes[ sCounter ] = particleIndex;
            ++sCounter;
          }

          //Particle in the box frame
          else if ( ( (halos[i].getX()-FOV) < partX ) && ( partX < (halos[i].getX()+FOV) ) &&
                    ( (halos[i].getY()-FOV) < partY ) && ( partY < (halos[i].getY()+FOV) ) &&
                    ( (halos[i].getZ()-FOV) < partZ ) && ( partZ < (halos[i].getZ()+FOV) ) ){
            boxIndexes[ bCounter ] = particleIndex;
            ++bCounter;
          }


          particleIndex = linkList[ particleIndex ];
        } //while
      }
      }   //Box loops
      }

      //Write the fits file for this halo
      writeFits( userInput ,
                  halos[i] ,
                 particles ,
                  N_sphere ,
                  N_box    ,
                  N_i1     ,
                  N_i2     ,
                  N_i3     ,
             sphereIndexes ,
                boxIndexes ,
                 i1Indexes ,
                 i2Indexes ,
                 i3Indexes );
//*/
exit(0);

    }     //If we found particles
  }       //Halo loop

}




void writeFits( inputInfo          userInput, //All the user info
                 haloInfo               halo, //The halo we are considering
                particlePosition particles[], //The full array of particles
                long long           N_sphere,
                long long           N_box   , //Number of particles in sphere set, box set, and integration length sets
                long long           N_i1    ,
                long long           N_i2    ,
                long long           N_i3    ,
                long long     *sphereIndexes, //Particle indexes in each set
                long long     *   boxIndexes,
                long long     *    i1Indexes,
                long long     *    i2Indexes,
                long long     *    i3Indexes){



  //Generate file names
  char temp[100];
  sprintf( temp,              "%sHalo%s_%li.FITS", (userInput.getParticleDir()).c_str(), (userInput.getCatType()).c_str(), halo.getID() );
  const std::string sphereFileName = temp;
  sprintf( temp,        "%sBox%s_%li_%04.1f.FITS", (userInput.getParticleDir()).c_str(), (userInput.getCatType()).c_str(), halo.getID(), userInput.getFOV() );
  const std::string    boxFileName = temp;

//  sprintf( temp, "%sBox%s_%li_%04.1f_i%06.1.FITS", (userInput.getParticleDir()).c_str(), (userInput.getCatType()).c_str(), halo.getID(), userInput.getFOV() );
  const std::string     i1FileName = temp;
//  sprintf( temp, "%sBox%s_%li_%04.1f_i%06.1.FITS", (userInput.getParticleDir()).c_str(), (userInput.getCatType()).c_str(), halo.getID(), userInput.getFOV() );
  const std::string     i2FileName = temp;
//  sprintf( temp, "%sBox%s_%li_%04.1f_i%06.1.FITS", (userInput.getParticleDir()).c_str(), (userInput.getCatType()).c_str(), halo.getID(), userInput.getFOV() );
  const std::string     i3FileName = temp;



//  sprintf( temp, "Box_%s_%s_%4.1f_i%5.1f.FITS", halo.getID(), userInput.getCatType(), userInput.getFOV() );
std::cout << sphereFileName << std::endl;
std::cout <<    boxFileName << std::endl;

std::cout<<"N's "<<
N_sphere<<" "<<
N_box   <<" "<<
N_i1    <<" "<<
N_i2    <<" "<<
N_i3    <<" "<<
std::endl;

printf("Halo: %li\nCenter: %6.2f %6.2f %6.2f\n", halo.getID(), halo.getX(), halo.getY(), halo.getZ());
for( int i = 0; i < N_sphere; ++i ){
  long long index = sphereIndexes[i];
  float         r = halo.getRm();
  float         x = particles[index].x_pos ;//- halo.getX();
  float         y = particles[index].y_pos ;//- halo.getY();
  float         z = particles[index].z_pos ;//- halo.getZ();
  printf("        %6.2f %6.2f %6.2f\n", x, y, z);
}
printf("Box:\n");
for( int i = 0; i < N_box; ++i ){
  long long index = boxIndexes[i];
  float         x = particles[index].x_pos - halo.getX();
  float         y = particles[index].y_pos - halo.getY();
  float         z = particles[index].z_pos - halo.getZ();
  printf("        %5.1f %5.1f %5.1f\n", x, y, z);
}



  //Number of elements in each direction, and total number of pixels
  int   N_pixels[2] = { userInput.getNPixlesH(), userInput.getNPixelsV() };
  long lN_pixels[2] = { N_pixels[0], N_pixels[1] };

  int   N_elements  = std::accumulate(&N_pixels[0],&N_pixels[2],1,std::multiplies<int>());


  //Our surface density array
  std::valarray<double> SD( N_elements );
  for ( int i = 0; i < N_elements; ++i ){
    SD[i] = 0;
  }

//find indexes based on x and y position indexes, don't even need z
//first do sphere, then box, then any integration lengths
//afterwards deallocate?
//then do LOS things

  float x_start = halo.getX() - userInput.getFOV() / 2.0;
  float y_start = halo.getY() - userInput.getFOV() / 2.0;

  float x_scale = N_pixels[0] / userInput.getFOV() ;
  float y_scale = N_pixels[1] / userInput.getFOV() ;


  for ( int i = 0; i < N_sphere; ++i ){

    long long index = sphereIndexes[i];

    float x_pos = particles[index].x_pos - x_start;
    float y_pos = particles[index].y_pos - y_start;

    int  hIndex = std::min( std::max(  int( x_pos * x_scale ) , 0 ), N_pixels[0] );
    int  vIndex = std::min( std::max(  int( y_pos * y_scale ) , 0 ), N_pixels[1] );

    SD[ vIndex * N_pixels[0] + hIndex ] += 1;
  }

//Consider particle mass

  //Fits stuff:

  {
    //Pointer to sphere Fits file
    std::auto_ptr<CCfits::FITS> pFits(0);

    //Create new FITS object
    pFits.reset( new CCfits::FITS( sphereFileName, DOUBLE_IMG, 2, lN_pixels) );

//    ( *pFits ).pHDU().addKey("OBJ",val,"DESC.");
    ( *pFits ).pHDU().write( 1, N_elements, SD);
  }


  //BOX, maybe we can make this a function?



  for ( int i = 0; i < N_sphere; ++i ){

    long long index =    boxIndexes[i];

    float x_pos = particles[index].x_pos - x_start;
    float y_pos = particles[index].y_pos - y_start;

    int  hIndex = std::min( std::max(  int( x_pos * x_scale ) , 0 ), N_pixels[0] );
    int  vIndex = std::min( std::max(  int( y_pos * y_scale ) , 0 ), N_pixels[1] );

    SD[ vIndex * N_pixels[0] + hIndex ] += 1;
  }

//Consider particle mass

  //Fits stuff:

  {
    //Pointer to sphere Fits file
    std::auto_ptr<CCfits::FITS> pFits(0);

    //Create new FITS object
    pFits.reset( new CCfits::FITS(    boxFileName, DOUBLE_IMG, 2, lN_pixels) );

//    ( *pFits ).pHDU().addKey("OBJ",val,"DESC.");
    ( *pFits ).pHDU().write( 1, N_elements, SD);
  }






/*
  //Pointer to the fits object
  std::auto_ptr<CCfits::FITS> pFits(0);

  try
  {
//    const std::string fileName("!atestfil.fit");
//    pFits.reset( new FITS(fileName , USHORT_IMG , naxis , naxes ) );
  }
  catch (FITS::CantCreate)
  {
    return -1;
  }

  //Change this?
  long fpixel(1);
  //  pFits->pHDU().addKey("EXPOSURE", exposure,     "Total Exposure Time");
  //  pFits->pHDU().write( fpixel, nelements, array);

*/


/*
  for (int i = 0; i < naxes[0]; ++i){
  for (int j = 0; j < naxes[1]; ++j){
    array[ i    + j * naxes[0] ] = 1000 * std::sin( ( float(i) / naxes[0] ) * M_PI / 2 ) * std::cos( ( float(j) / naxes[1] ) * M_PI / 2 ) ;
  }
  }


here
i horizontal, j vertical
both start bottom right in gimp viewer
*/

}
