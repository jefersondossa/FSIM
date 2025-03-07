#ifndef ASSEMBLE_H
#define ASSEMBLE_H

#include "Analysis.h"
#include "LinearAnalysis.h"

namespace Assemble{
    void Monomodel(Analysis *fAnalysis, int mesh = 0, int64_t startDOF = 0);
    void MonomodelMatrix(Analysis *fAnalysis, int mesh = 0, int64_t startDOF = 0);
    void MonomodelVector(Analysis *fAnalysis, int mesh = 0, int64_t startDOF = 0);
    void Arlequin(Analysis *fAnalysis);
    void ArlequinMatrix(Analysis *fAnalysis);
    void ArlequinVector(Analysis *fAnalysis);
    void Coupling(Analysis *fAnalysis, int64_t startDOF = 0);
    void CouplingMatrix(Analysis *fAnalysis, int64_t startDOF = 0);
    void CouplingVector(Analysis *fAnalysis, int64_t startDOF = 0);
    void LagrangeMultiplierDOF(Analysis *fAnalysis, int mesh = 0, int64_t startDOF = 0);
    void LagrangeMultiplierDOFMatrix(Analysis *fAnalysis, int mesh = 0, int64_t startDOF = 0);
    void LagrangeMultiplierDOFVector(Analysis *fAnalysis, int mesh = 0, int64_t startDOF = 0);
    void stabilizeArlequin(Analysis *fAnalysis, std::vector<MatrixDouble> &Stiffness, std::vector<VecDouble> &Rhs, int64_t &element);
};




#endif