#include <iostream>
#include <cstring>

#ifndef HE_CLASSES
#define HE_CLASSES

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

#endif // HE_CLASSES
