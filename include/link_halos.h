

#ifndef LINK_HALOS
#define LINK_HALOS


//Make the link list of nearby particles
void makeLinkList( particlePosition particle );

//Sets the min/max values of the particle positions
void setMinMaxParticles( particlePosition particle[],  //Particles to find min/max of
                         inputInfo       &userInput ); //Number of particles in the array

#endif // LINK_HALOS
