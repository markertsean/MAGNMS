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

void makeLinkList( particlePosition particle ){


}
