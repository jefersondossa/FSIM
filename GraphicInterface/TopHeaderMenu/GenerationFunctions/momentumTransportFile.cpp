#include "WindowConstructor.h"

void WindowConstructor::fmomentumTransportFile(){

    fstream mT_file;

    //Lets open the Ufile and write his header

    string mTheader;
    fstream headerFrommT;

    headerFrommT.open("../GraphicInterface/TopHeaderMenu/GenerationFunctions/headerfrommomentumTransport.txt",ios::in);
    if(headerFrommT.is_open()){ //Copies text from headerFrommT to a string mTheader.
      while(getline(headerFrommT, Line)){
        mTheader += Line + "\n";
      }
      headerFrommT.close();
    }

    mT_file.open((constantFiles[0]),ios::app);
    if (mT_file.is_open()){

        mT_file << mTheader;
        mT_file.close();
    }
}