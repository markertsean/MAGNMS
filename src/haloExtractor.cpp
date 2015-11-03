#include <iostream>
#include <cstring>


//Read halo catalog
//Select halos within criterion
//Write useable catalog, if one doesn't exist
//Create directories for infofile, particle positions
//Go through particle catalog, extracting

//Need read functions for particle file, halo catalog
//Identify catalog/particle file type

//Function to to create directories using file names

class inputInfo{

    void setInputCatalog( std::string inpS ) { inputCatalog   = inpS; }
    void setInputPart   ( std::string inpS ) { inputPartFiles = inpS; }
    void setparticleDir ( std::string inpS ) { particleDir    = inpS; }
    void setheaderDir   ( std::string inpS ) {   headerDir    = inpS; }

    std::string getInputCatalog() { return inputCatalog   ; }
    std::string setInputPart   () { return inputPartFiles ; }
    std::string setparticleDir () { return particleDir    ; }
    std::string setheaderDir   () { return   headerDir    ; }

  private:
    std::string inputCatalog   = " ";
    std::string inputPartFiles = " ";

    std::string particleDir    = " ";
    std::string   headerDir    = " ";
};


int main( int arg, char ** argv ){

//Need to read halo catalog, particle catalog name, optionally directory name
//Read the halo catalog

  inputInfo userInput;


}
