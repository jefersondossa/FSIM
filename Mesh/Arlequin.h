//------------------------------------------------------------------------------
// 
//                   Jeferson W D Fernandes and Rodolfo A K Sanches
//                             University of Sao Paulo
//                           (C) 2017 All Rights Reserved
//
// <LicenseText>
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//----------------------------------ARLEQUIN------------------------------------
//------------------------------------------------------------------------------

#ifndef ARLEQUIN_H
#define ARLEQUIN_H

#include "IntegrationQuadrature11.h"
#include <map>
#include <set>
#include "Analysis.h"
#include "CompMesh.h"
#include "Element.h"
// #include "VTUGenerator.h"


class Analysis;
class Element;
class CompMesh;
// class VTUGenerator;

enum ArlequinStabType{ENoStab, EOption1, EOption2, EOption3, EOption4, EOption5};

/// Mounts the overlapping mesh problem for solving the incompressible flow problem
class Arlequin{
public:
    std::set<int64_t> fGluingElementIndex;

    std::map<int64_t,int64_t> fNodeLocalToElementGlobal;
    std::map<int64_t,VecDouble> fNodeLocalToXsiGlobal;

    std::map<int64_t,VecDouble> fLocalIntPointToGlobalElement;
    std::map<int64_t,MatrixDouble> fLocalIntPointToGlobalXsi;

    std::map<int64_t,std::set<int64_t>> fGlobalElToLocalEl;
    std::set<int> fGlueMatID;
private:
    double fGlueZoneThickness = 0.0;
    double fArlequinEpsilon = 1.e-3;

    int numElemGlueZoneFine;
    int numElemGlueZoneCoarse;
    int numNodesGlueZoneFine;
    int numNodesGlueZoneCoarse;
    bool fFirstSearch = true;
    bool fInvertSignaledDistance = false;

    /// @brief Data structure relating signaled distance of local model nodes
    std::map<int64_t,double> fLocalSignaledDistance;
    /// @brief Data structure relating signaled distance of global model nodes
    std::map<int64_t,double> fGlobalSignaledDistance;
    
    /// @brief Arlequin stabilization type
    ArlequinStabType fArlequinStab = ArlequinStabType::ENoStab;
    std::vector<CompMesh *> fMeshVector;

    double fK0;
    double fK1;
public:
    /// @brief Default constructor
    Arlequin() = default;

    /// @brief Default destructor
    ~Arlequin() = default;

    /// @brief Arlequin model Constructor
    /// @param meshvec Mesh vector
    /// @param k0 L2 coupling
    /// @param k1 H1 coupling
    /// @param stab Arlequin stabilization option
    Arlequin(std::vector<CompMesh *> &meshvec, double k0, double k1, ArlequinStabType stab = ArlequinStabType::ENoStab);

    /// @brief Sets the glue material ids
    /// @param glue glue material id
    void SetGlueIds(std::set<int> &glue){
        fGlueMatID = glue;
    }

    /// @brief Returns the mesh vector
    /// @return mesh vector
    std::vector<CompMesh *> &MeshVec(){return fMeshVector;}

    /// @brief Sets the Arlequin stabilization option
    /// @param stab stabilization type
    void setArlequinStabilization(ArlequinStabType stab){
        fArlequinStab = stab;
    }

    /// @brief returns the Arlequin stabilization option 
    /// @return Arlequin stabilization option
    ArlequinStabType &getArlequinStabilization(){return fArlequinStab;}

    /// @brief Compute and store the element boxes for improving the correspondence searching process
    void SetElementBoxes();

    /// @brief Defines the gluing (or coupling) zone
    void setCouplingZone();

    /// @brief Compute and store the signaled distance function
    void setSignaledDistance();

    /// @brief Sets the energy weight function 
    /// @param double reference value for computing the energy weight function
    void setWeightFunction(double val);

    /// @brief Sets the nodal correspondence of the fine to the coarse models
    void setNodalCorrespondenceFine();

    /// @brief Searchs the node correspondence in a computational mesh
    /// @param x point
    /// @param cmesh the computational mesh
    /// @param elCorr correpondence element
    /// @param xsiCorr correspondence Xsi
    /// @param elSearch the first guess element
    void searchNodeCorrespondence(VecDouble &x, CompMesh *cmesh,
                                  int &elCorr, VecDouble &xsiCorr, int elSearch);

    void CreateGlobalCouplingElements();

    void SetUp(){
        SetElementBoxes();
        setSignaledDistance();

        //Construct the glue zone based on some defined criterion
        setCouplingZone();
        //Computes the Weight function for all the finite elements
        setWeightFunction(16.);
        
        setNodalCorrespondenceFine();
        CreateGlobalCouplingElements();
    };

    void SetGlueZoneThichkess(double val){
        fGlueZoneThickness = val;
    }
    void SetEpsilon(double val){
        fArlequinEpsilon = val;
    }

    double GlobalWeightFunction(double dist);
    double LocalWeightFunction(double dist);

    void ComputeTriangleIntersections(int64_t iEl, std::set<int64_t> &elIntersected);

    double LocalNodeSignaledDistance(int64_t node){
        return fLocalSignaledDistance[node];
    }
    double GlobalNodeSignaledDistance(int64_t node){
        return fGlobalSignaledDistance[node];
    }

    double ShortestDistance(VecDouble &point, VecDouble &startSeg, VecDouble &endSeg);

    void InvertSignaledDistance(){
        fInvertSignaledDistance = true;
    }
};


#endif
