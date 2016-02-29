
#include "halo_extraction_classes.h"

#ifndef LINK_HALOS
#define LINK_HALOS



//Sets the min/max values of the particle positions
void setMinMaxParticles( particlePosition particle[],  //Particles to find min/max of
                         inputInfo       &userInput ); //Number of particles in the array

//Locates halos that could be in our box with a full FOV
unsigned long findBoxHalos( inputInfo &userInfo   ,
                             haloInfo  allHalos[] ,
                             haloInfo  boxHalos[] ,
                                short  runNum     );

//Make the link list of nearby particles
void makeLinkList( inputInfo        userInfo   ,   //Contains the global info needed
                   particlePosition particle[] ,   //Array of the particles
                   long long        myList  [] ,   //List pointing to next neighbor particle
                   long long        myLabel [] );  //Label points to the last particle in list for index


void linkHaloParticles( inputInfo userInput   ,
                         haloInfo     halos[] ,
                 particlePosition particles[] ,
                       long long  labelList[] ,
                       long long   linkList[] );


void writeImage( inputInfo          userInput  , // All the user info
                  haloInfo               halo  , // The halo we are considering
                 particlePosition particles[]  , // The full array of particles
                 long long          N_sphere   ,
                 long long          N_box      , // Number of particles in sphere set, box set, and integration length sets
                 long long          N_integ[]  ,
                 long long       maxN_integ    , // Maximum number on particles in any integIndexes bin
                 long long          N_integBins, // Number of bins (steps) of integration
                 double          integLengths[], // Contains the integration lengths
                 long long     *sphereIndexes  , // Particle indexes in each set
                 long long     *   boxIndexes  ,
                 long long     * integIndexes  );


int  writeFits( const std::string     fileName    ,  // File name to write
                int                   N_pixels[]  ,  // N_pixels in each direction
                int                   N_pixelsTot ,  // Total number of pixels
                std::valarray<double>        SD   ,  // SD array to use, passed because we will keep adding to the array
                long long             N_indexes   ,  // Number of indexes in set to add
                long long               indexes[] ,  // Indexes of the set
                particlePosition      particles[] ,  // All the particles
                haloInfo                   halo   ,  // Central halo
                inputInfo             userInput   ); // All the user information


#endif // LINK_HALOS
