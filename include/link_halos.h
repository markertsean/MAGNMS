
#include "halo_extraction_classes.h"

#ifndef LINK_HALOS
#define LINK_HALOS


//Make the link list of nearby particles
void makeLinkList( particlePosition particle );

//Sets the min/max values of the particle positions
void setMinMaxParticles( particlePosition particle[],  //Particles to find min/max of
                         inputInfo       &userInput ); //Number of particles in the array

//Locates halos that could be in our box with a full FOV
unsigned long findBoxHalos( inputInfo &userInfo   ,
                             haloInfo  allHalos[] ,
                             haloInfo  boxHalos[] ,
                                short  runNum     );

#endif // LINK_HALOS
