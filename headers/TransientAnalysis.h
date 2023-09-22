#ifndef TRANSIENTANALYSIS_H
#define TRANSIENTANALYSIS_H

#include "LinearAnalysis.h"
#include "NonLinearAnalysis.h"

class TransientAnalysis : public LinearAnalysis, public NonLinearAnalysis
{
private:
    /* data */
public:
    // TransientAnalysis(/* args */);
    // ~TransientAnalysis();
};


// TransientAnalysis::~TransientAnalysis()
// {
// }


#endif