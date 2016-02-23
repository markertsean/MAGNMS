#include <iostream>
#include <fstream>
#include <cstring>
#include <sys/stat.h>

#include <CCfits/CCfits>
#include <cmath>

#include "halo_extraction_classes.h"
#include "read_files.h"
#include "link_halos.h"

//Read halo catalog
//Select halos within criterion
//Write useable catalog, if one doesn't exist
//Create directories for infofile, particle positions
//Go through particle catalog, extracting

//Need read functions for particle file, halo catalog
//Identify catalog/particle file type

//Function to to create directories using file names




int main( int arg, char ** argv ){

  //////////////////////////////////////
  //////////User info read in///////////
  //////////////////////////////////////


  //Stores the user input from the input file
  inputInfo userInput;

  //User can specify input file via command line, default: extractInfo.dat
  if ( arg > 1 ){
    userInput.setReadFile( argv[1] );
  }
  printf("\n Using input file: %s\n", (userInput.getReadFile()).c_str());

  //Attempts to read the input file, and displays the resultant info
  if ( readUserInput( userInput.getReadFile(), userInput ) ){
    printf("  Halo Catalog       : %s\n  Particle File      : %s\n  Header Directory   : %s\n  Particle Directory : %s\n\n",
    (userInput.getInputCatalog()).c_str(),
    (userInput.getInputPart   ()).c_str(),
    (userInput.getHeaderDir   ()).c_str(),
    (userInput.getParticleDir ()).c_str());
  }
  else {
    printf("Error in reading input file, aborting\n\n");
    exit(1);
  }


  ///////////////////////////////////////////
  ////////Read in the catalog file///////////
  ///////////////////////////////////////////



  //Read the number of valid halos, for allocation
  unsigned long N_halos =0;
  {
    haloInfo myHalos[2];
    N_halos = readCatalog( myHalos, userInput, N_halos );
    printf(" Number of halos read: %lu\n",N_halos);
  }

//  haloInfo myHalos[N_halos];
  haloInfo *myHalos = new haloInfo[N_halos];

  {
    unsigned long old_N_halos = N_halos;
    N_halos = readCatalog( myHalos, userInput, N_halos );
    if ( old_N_halos != N_halos ){
      printf("\nError: Read mismatch\nN_halos allocated: %lu\nN_halos read in: %lu\n\n", old_N_halos, N_halos);
      exit(1);
    }
  }
  userInput.setNumHalos( N_halos );


  printf("\n Catalog read in complete\n\n");


  /////////////////////////////////////////
  //////Write header/part directory////////
  /////////////////////////////////////////

  {
    struct stat sb;
    int maxStrLength=400;
    char str[maxStrLength];

    //If header directory exists, we use it
    //If not, create it
    sprintf( str, "mkdir %s", (userInput.getHeaderDir()).c_str() );
    if ( stat((userInput.getHeaderDir()).c_str(), &sb) == 0 ){
      std::cout << " Found directory: " << userInput.getHeaderDir() << std::endl;;
    } else {
      std::cout << " Writing directory: " << userInput.getHeaderDir() << std::endl;
      system(str);
    }

    sprintf( str, "mkdir %s", (userInput.getParticleDir()).c_str() );
    if ( stat((userInput.getParticleDir()).c_str(), &sb) == 0 ){
      std::cout << " Found directory: " << userInput.getParticleDir() << std::endl;;
    } else {
      std::cout << " Writing directory: " << userInput.getParticleDir() << std::endl;
      system(str);
    }

    printf("\n\n");
  }



  /////////////////////////////////////////
  ///////////Read in particles/////////////
  /////////////////////////////////////////

  //setPartFile will locate the particle files
  // if the file exists, or create the file
  // based on the PMSS file it thinks we should use
  //If the read in fails, will return a 0 for num particles

  long long numParticles = 0;

  numParticles = setPartFile( userInput );

  printf("\n");
  if ( numParticles == 0 ){

    printf(" Error reading particle file:%s\n", (userInput.getInputPart()).c_str() );
    printf("Specify the file using the variable \"partFile\", using the variable \"snapNum\", and making sure you are in the proper directory\n\n");
    exit(1);

  }
  userInput.setNumParticles( numParticles );
  //Once particles are confirmed, read into this C side
  //particlePosition particle[numParticles];                       //This method segfaults due to compiler limit

  particlePosition *particle = new particlePosition[numParticles]; //This does not segfault


  readParticle( userInput, particle );


  /////////////////////////////////////////
  ////////Generate link lists//////////////
  /////////////////////////////////////////

  std::cout << " Setting up cells..." << std::endl;

  //Determines range of diff coordinates
  setMinMaxParticles( particle, userInput );

  //Find the halos who's FOV is entirely inside the box
  {
  //Dummy array, useless
  haloInfo duhalos[2];

  N_halos = findBoxHalos( userInput, myHalos, duhalos, 0 ); //Returns the number of halos in the box
  }

  haloInfo halos[N_halos];                                  //Allocates our new halo array
  N_halos = findBoxHalos( userInput, myHalos,   halos, 1 ); //Fills in our new halo array
  delete [] myHalos;                                        //Deletes old halo array


  if ( !userInput.setBox( userInput.getFOV() / 2.0 ) ){
    std::cout << " Error setting boxes "  << std::endl;
    exit(1);
  }

  std::cout << " Done."     << std::endl  << std::endl;


  std::cout << " Generating link list..." << std::endl;

  //Make link list between particles
  long long  *linkList = new long long [ userInput.getNumParticles() ];
  long long *labelList = new long long [ userInput.getNtotCell()     ];

  makeLinkList( userInput, particle, linkList, labelList );
  std::cout << " Done."     << std::endl  << std::endl;


  /////////////////////////////////////////
  ////////Link halos, write files//////////
  /////////////////////////////////////////


//Loop over halos, find boxes to check with
//Z axis should check against long integration lengths
//Loop once to find, generate link lists for each set
//sphere, box, integration lengths (if applicable)
//Integration length images first, write file
//Box, write file
//Rotate sphere for different views, write file


//For now do arrays, figure out the fits writing later

//32 boxes
  std::cout << " Generating SD boxes..." << std::endl;
  linkHaloParticles( userInput, halos, particle, labelList, linkList );
  std::cout << " Done."     << std::endl  << std::endl;



exit(0);

using namespace CCfits;
{

  std::auto_ptr<FITS> pFits(0);
  long   naxis    = 2;
  long   naxes[2] = { 300, 200 };
  long  nelements = std::accumulate(&naxes[0],&naxes[naxis],1,std::multiplies<long>());
  std::valarray<int> array(nelements);

  for (int i = 0; i < naxes[0]; ++i){
  for (int j = 0; j < naxes[1]; ++j){
    array[ i    + j * naxes[0] ] = 1000 * std::sin( ( float(i) / naxes[0] ) * M_PI / 2 ) * std::cos( ( float(j) / naxes[1] ) * M_PI / 2 ) ;
  }
  }

/*
here
i horizontal, j vertical
both start bottom right in gimp viewer
*/

  try
  {
    const std::string fileName("!atestfil.fit");
    pFits.reset( new FITS(fileName , USHORT_IMG , naxis , naxes ) );
  }
  catch (FITS::CantCreate){  return -1; }


  long fpixel(1);
  long exposure(1500);
  std::complex<float> omega(std::cos(2*M_PI/3.),std::sin(2*M_PI/3));


  pFits->pHDU().addKey("EXPOSURE", exposure,     "Total Exposure Time");
  pFits->pHDU().addKey("OMEGA"   ,    omega," Complex cube root of 1 ");
  pFits->pHDU().write( fpixel, nelements, array);
}

exit(0);
/*
// Create a FITS primary array containing a 2-D image declare axis arrays.
long naxis = 2;
long naxes[2] = { 300, 200 };

// references for clarity.
long& vectorLength = naxes[0];
long& numberOfRows = naxes[1];
long  nelements              ;
// Find the total size of the array.
// this is a little fancier than necessary ( It’s only calculating naxes[0]*naxes[1]) but it demonstrates use of the C++ standard library accumulate algorithm.
nelements = std::accumulate(&naxes[0],&naxes[naxis],1,std::multiplies<long>());

//Gen values
// create a dummy row with a ramp. Create an array and copy the row to
// row-sized slices. [also demonstrates the use of valarray slices].
// also demonstrate implicit type conversion when writing to the image:
// input array will be of type float.
std::valarray<int> row(vectorLength);
std::valarray<int> array(nelements);
for (long j = 0; j < vectorLength; ++j) row[j] = j; //Number of columns
for (int  i = 0; i < numberOfRows; ++i)
{
array[
std::slice(vectorLength*static_cast<int>(i),vectorLength,1)
] = row + i;
}

// Creates 300x300 array
std::vector<long> extAx(2,300);

// create some data for the image extension.
                                  //     range to do, initial val, multiply them, just the number of elements
long extElements = std::accumulate( extAx.begin(), extAx.end(), 1, std::multiplies<long>());
//Gen values
// Array of floats
std::valarray<float> ranData(extElements);
const float PIBY (M_PI/150.);
//Generate cos values across rows
for ( int jj = 0 ; jj < extElements ; ++jj)
{
float arg = PIBY*jj;
ranData[jj] = std::cos(arg);
}

//add two keys to the primary header, one long, one complex.
long exposure(1500);
std::complex<float> omega(std::cos(2*M_PI/3.),std::sin(2*M_PI/3));



// declare auto-pointer to FITS at function scope. Ensures no resources leaked if something fails in dynamic allocation.
std::auto_ptr<FITS> pFits(0);

// Create new image
try
{
  // overwrite existing file if the file already exists.
  const std::string fileName("!atestfil.fit");
  // Create a new FITS object, specifying the data type and axes for the primary image. Simultaneously create the corresponding file.
  // this image is unsigned short data, demonstrating the cfitsio extension to the FITS standard.
  pFits.reset( new FITS(fileName , USHORT_IMG , naxis , naxes ) );
}
catch (FITS::CantCreate)
{
  // ... or not, as the case may be.
  return -1;
}

// Create image header?
std::string newName ("NEW-EXTENSION");
ExtHDU* imageExt = pFits->addImage(newName,FLOAT_IMG,extAx);

long fpixel(1);


//Array is [300,200], ranData is [300,300] extAxis, number 2?

//pfits the pointer to the fits thing
//ExtHDU * imageExt = pFIts->addImage(<name?>,<data type?>,<array?>)
//array and ranData hold data
//wrote ranData to imageExt
//wrote stuff to header
//wrote data to image


//std::auto_ptr<FITS> pFits(0);
//try
//{
//  const std::string fileName("!atestfil.fit");
//  pFits.reset( new FITS(fileName , USHORT_IMG , naxis , naxes ) );
//}
//catch (FITS::CantCreate){  return -1; }
//long fpixel(1);
//pFits->pHDU().addKey("EXPOSURE", exposure,     "Total Exposure Time");
//pFits->pHDU().addKey("OMEGA"   ,    omega," Complex cube root of 1 ");
//pFits->pHDU().write( fpixel, nelements, array);



// write the image extension data: also demonstrates switching between
// HDUs.
imageExt->write(fpixel,extElements,ranData);

pFits->pHDU().addKey("EXPOSURE", exposure,     "Total Exposure Time");
pFits->pHDU().addKey("OMEGA"   ,    omega," Complex cube root of 1 ");

// The function PHDU& FITS::pHDU() returns a reference to the object representing
// the primary HDU; PHDU::write( <args> ) is then used to write the data.
pFits->pHDU().write( fpixel, nelements, array);

// PHDU’s friend ostream operator. Doesn’t print the entire array, just the
// required & user keywords, and is provided largely for testing purposes
// [see readImage() for an example of how to output the image array to a stream].
std::cout << pFits->pHDU() << std::endl;
return 0;
*/



  //Catalog
  //ID number
  //Halo position
  //Particle number from catalog
  //Mass of particles
  //Halo mass
  //Rmax
  //Concentration
  //b/a ratio
  //c/a ratio
  //
  //Rotation vector for along LOS
  //Rotation vector for perp to LOS
  //Number particles using
  //Number particles in box
  //Whether could create LOS cones
  //Largest LOS cones

  //Have input files and directories to write to.
  //Need to:
  //        Write overarching header file, with used halos
  //        Read in particles, maybe use link list to speed up process?



}
