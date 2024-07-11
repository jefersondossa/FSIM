#include "WindowConstructor.h"

namespace fs = std::filesystem;

void WindowConstructor::fFoamFile(){

    fOFfolderGenerator();
    fUFile();
    fpFile();
    fphysicalPropertiesFile();
    fmomentumTransportFile();
}