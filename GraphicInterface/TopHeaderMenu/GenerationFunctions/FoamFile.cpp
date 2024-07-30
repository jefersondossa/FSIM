#include "WindowConstructor.h"

namespace fs = std::filesystem;

void WindowConstructor::fFoamFile(){

    fOFfolderGenerator();
    fUFile();
    fpFile();
    fkFile();
    fEpsilonFile();
    fnutFile();
    fphysicalPropertiesFile();
    fmomentumTransportFile();
    fcontrolDictFile();
    ffv();
    fRunFoam();
}

