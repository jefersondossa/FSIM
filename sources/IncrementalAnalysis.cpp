#include "IncrementalAnalysis.h"

void IncrementalAnalysis::Run(){

    int iStep = 0;

    while (iStep < fNSteps)
    {   
        NonLinearAnalysis::Run();
    }
}
