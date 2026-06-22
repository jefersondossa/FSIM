#ifndef MEFGGlobalLocalTools_H
#define MEFGGlobalLocalTools_H

#include "Element.h"
#include "CompMesh.h"
#include "GlobalLocalEnrichment.h"
#include "MixedGlobalLocalEnrichment.h"

class CompMesh;
class Element;

namespace MEFGGlobalLocalTools{

    void LocalToGlobalCorrespondence(CompMesh *cmeshG, CompMesh *cmeshL, 
                                     std::map<int,int> &globalElementCorrespondence, 
                                     std::map<int, MatrixDouble> &globalNodeCorrespondence,
                                     int overlappingRegion, int overlappingNHDirichletBoundary = -1);
    
    
    void LocalToGlobalCorrespondenceBoundary(CompMesh *cmeshG, CompMesh *cmeshL,
                                           std::map<int,int> &globalElementCorrespondence,
                                           std::map<int, MatrixDouble> &globalNodeCorrespondence,
                                           int overlappingNHNeumannBoundary = -1);

    void CreateEnrichedModel(CompMesh *cmeshG, CompMesh *cmeshL, 
                            std::map<int,int> &globalElementCorrespondence, 
                            std::map<int, MatrixDouble> &globalNodeCorrespondence,
                            std::map<int,int> &enrichedConnects,
                            int overlappingNHNeumannBoundary,
                            GlobalLocalEnrichment *globalLocal);

    void CreateMixedEnrichedModel(CompMesh *cmeshG, CompMesh *cmeshL, 
                                std::map<int,int> &globalElementCorrespondence, 
                                std::map<int, MatrixDouble> &globalNodeCorrespondence,
                                std::map<int,int> &enrichedConnects,
                                int overlappingNHNeumannBoundary,
                                MixedGlobalLocalEnrichment *globalLocal);
};


#endif