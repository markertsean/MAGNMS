#include <iostream>
#include <cstring>
#include "halo_extraction_classes.h"


//Reads input files and optional directories to write files to (otherwise default)
bool readUserInput( std::string fileName, inputInfo &myInput ){

  FILE *pFile;
  char   inpC1[50],inpC2[50];
  pFile = fopen( fileName.c_str() ,"r");

  //If we can open the file, cycle through lines
  if ( pFile != NULL ){
    while ( fscanf(pFile,"%s%s",inpC1,inpC2) != EOF ){

      std::string inpS = std::string(inpC2);

      if      ( strcmp( inpC1, "haloFile" ) == 0 ){
        myInput.setInputCatalog( inpS  );
      }
      else if ( strcmp( inpC1, "partFile" ) == 0 ){
        myInput.setInputPart   ( inpS  );
      }
      else if ( strcmp( inpC1, "headDir"  ) == 0 ){
        myInput.setParticleDir ( inpS  );
      }
      else if ( strcmp( inpC1, "partDir"  ) == 0 ){
        myInput.setHeaderDir   ( inpS  );
      }

    }
  }
  //If file doesn't open, throw error
  else {
    std::cout << "Could not open input file: " << fileName << std::endl;
    return false;
  }

  fclose(pFile);

  if ( myInput.getInputCatalog() == "" || myInput.getInputPart() == "" ){
    printf("Missing required input files:\n   Halo catalog (haloFile) = %s\n   Particle file (partFile) = %s\n\n",
          (myInput.getInputCatalog()).c_str(),(myInput.getInputPart()).c_str());
    return false;
  }

  myInput.setDirectory();

  return true;
}


bool readCatalog(){
//Pass user info
//Pass halo array
/*
bigMD2.5Om27
 A    = 1.00000
 I    =    0 Nrow =    0 Ngrid=    0
 Omega_0=  0.30711 Omega_L=  0.69288 hubble =  0.67770 buffer width (Mpch) =  5.00000
 Force softening (kpch)               =   5.000
 Minimum distance between halos (kpch)=  25.000
 Number of neighbours for  particles  = 20
 Number of buffered particles         =  20000
 Split factor for adaptive mesh       =  8
 Number of radial bins                =   30
 Step in radius dlog10                = 0.05000
 Mass of smallest particle (Msunh)    =  1.5054E+09
 Subhalos rejection: virial ratio     =   1.200
 offset parameter limit               =   0.200
 radius overshoot factor              =   1.500
 Overdensity limit                    =  651.22
     XYZ(Mpch)                          Vxyz(km/s)                  Mbound     Mtot/Msunh    Rvir(kpch) Vrms(km/s) Vcirc(km/s)       Nhalo  Cvir    Nparticles  Distinct/Sub   Xoff  2K/Ep-1   Lambda   RadRMS/kpch
  b/a  c/a MajorAxis:  x      y      z
     1.5268     3.2015     1.0790       -52.14    284.31    181.71  3.7635E+10  3.7635E+10  54.730       52.04       59.22               1   5.872       25.00               0  0.1328      8.0193E-02  5.4501E-02
 30.41      0.6417      0.6194     -0.3962      0.1408      0.9073
*/
  return true;
}
