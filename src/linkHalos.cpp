#include <math.h>


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
                   unsigned long    myList  [] ,   //List pointing to next neighbor particle
                   unsigned long    myLabel [] ){  //Label points to the last particle in list for index

  double cell = 25;

//Can we leave it as is, or need to account for offsets?

  int Nlx = (                      userInfo.getXmin() ) / cell ; //Minimum box numbers
  int Nly = (                      userInfo.getYmin() ) / cell ;
  int Nlz = (                      userInfo.getZmin() ) / cell ;
  int Nrx = ( userInfo.getXmax() - userInfo.getXmin() ) / cell ; //Maximum box numbers
  int Nry = ( userInfo.getYmax() - userInfo.getYmin() ) / cell ;
  int Nrz = ( userInfo.getZmax() - userInfo.getZmin() ) / cell ;


  //Fortran has -1, check if this matters, we are using unsigned long however
  for ( unsigned long i = 0 ; i < userInfo.getNumParticles() ; ++i ){
    myList[i] = 0;
  }

  for ( int i = 0; i < Nrx ; ++i ){  //x
  for ( int j = 0; j < Nry ; ++j ){  //y
  for ( int k = 0; k < Nrz ; ++k ){  //z
    myLabel[ i + j * Nrx + k * Nrx * Nry ] = 0;
  }
  }
  }

  //Loop over each particle, generating the link list
  for ( unsigned long i = 0; i < userInfo.getNumParticles() ; ++i ){

    //Find which box the particle is in, make sure we stay in the box
    int xIndex = std::min( std::max( Nlx, int( ( particle[i].x_pos - userInfo.getXmin() ) / cell ) ), Nrx );
    int yIndex = std::min( std::max( Nly, int( ( particle[i].y_pos - userInfo.getYmin() ) / cell ) ), Nry );
    int zIndex = std::min( std::max( Nlz, int( ( particle[i].z_pos - userInfo.getZmin() ) / cell ) ), Nrz );


    unsigned long lI = xIndex + yIndex * Nrx + zIndex * Nrx * Nry; //Label index


    myList[  i  ] = myLabel[ lI ]; //List previous particle in the same box
    myLabel[ lI ] = i;             //label stores last particle recorded in the box
  }

}

/*
!--------------------------------------------------------------
!               Make linker lists of particles in each cell
      SUBROUTINE List
!--------------------------------------------------------------
      Cell = rmax
         Nmx         = xl/Cell - 1
         Nmy         = yl/Cell - 1
         Nmz         = zl/Cell - 1
         Nbx         = xr/Cell + 1
         Nby         = yr/Cell + 1
         Nbz         = zr/Cell + 1
         Allocate(Lst(Ntot))
         Allocate(Label(Nmx:Nbx,Nmy:Nby,Nmz:Nbz))
         CALL OMP_SET_NUM_THREADS(num_threads)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (i)
             Do i=1,Ntot
                Lst(i)=-1
             EndDo
         CALL OMP_SET_NUM_THREADS(num_threads)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (i,j,k)
             Do k=Nmz,Nbz
             Do j=Nmy,Nby
             Do i=Nmx,Nbx
                Label(i,j,k)=0
             EndDo
             EndDo
             EndDo
      Do jp=1,Ntot
         i=Ceiling(Xp(jp)/Cell)-1
         j=Ceiling(Yp(jp)/Cell)-1
         k=Ceiling(Zp(jp)/Cell)-1
         i=MIN(MAX(Nmx,i),Nbx)
         j=MIN(MAX(Nmy,j),Nby)
         k=MIN(MAX(Nmz,k),Nbz)
         Lst(jp)      =Label(i,j,k)
         Label(i,j,k) =jp
      EndDo
    end SUBROUTINE List
*/
