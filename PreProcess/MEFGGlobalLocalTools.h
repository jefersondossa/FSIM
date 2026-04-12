#ifndef MEFGGlobalLocalTools_H
#define MEFGGlobalLocalTools_H

#include "Element.h"
#include "CompMesh.h"
#include "GlobalLocalEnrichment.h"

class CompMesh;
class Element;

namespace MEFGGlobalLocalTools{

    void LocalToGlobalCorrespondence(CompMesh *cmeshG, CompMesh *cmeshL, 
                                     std::map<int64_t,int64_t> &globalElementCorrespondence, 
                                     std::map<int64_t, MatrixDouble> &globalNodeCorrespondence,
                                     int overlappingRegion, int overlappingNHDirichletBoundary = -1);
    
    
    void LocalToGlobalCorrespondenceBoundary(CompMesh *cmeshG, CompMesh *cmeshL,
                                           std::map<int64_t,int64_t> &globalElementCorrespondence,
                                           std::map<int64_t, MatrixDouble> &globalNodeCorrespondence,
                                           int overlappingNHNeumannBoundary = -1);

    void CreateEnrichedModel(CompMesh *cmeshG, CompMesh *cmeshL, 
                            std::map<int64_t,int64_t> &globalElementCorrespondence, 
                            std::map<int64_t, MatrixDouble> &globalNodeCorrespondence,
                            std::map<int64_t,int64_t> &enrichedConnects,
                            int overlappingNHNeumannBoundary,
                            GlobalLocalEnrichment *globalLocal);

    void CreateMixedEnrichedModel(CompMesh *cmeshG, CompMesh *cmeshL, 
                                std::map<int64_t,int64_t> &globalElementCorrespondence, 
                                std::map<int64_t, MatrixDouble> &globalNodeCorrespondence,
                                std::map<int64_t,int64_t> &enrichedConnects,
                                int overlappingNHNeumannBoundary,
                                GlobalLocalEnrichment *globalLocal);
};


#endif