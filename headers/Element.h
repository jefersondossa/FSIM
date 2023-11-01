#ifndef ELEMENT
#define ELEMENT

#include "DataTypes.h"
#include "CompMesh.h"
#include "IntPointData.h"

class CompMesh;

class Element{
protected:
    CompMesh *fMesh;
    VecInt        fConnect; //Velocity mesh connectivity 
    int64_t       fIndex;             //Element index

    int           fSideInBoundary;

    VecDouble intPointWeightFunction;
    VecDouble fIntPointDistFunction;
    VecDouble intPointWeightFunctionPrev;

    MatrixDouble fIntPointCoordinates;
    int DEG;

    VecDouble     xK, XK;
    
    

    std::vector<int64_t> fNeighborElements;

public:

    bool          FSIInterface;  

public:
    Element() = default;

    Element(int64_t index, VecInt &connect, CompMesh* mesh){
        
        // fMesh = mesh;
        // fConnect.resize(fMesh->NElNodes());
        // fIndex = index;
        // for (int i = fMesh->NElNodes(); i--; ) fConnect[i] = connect[i];
        // DEG = fMesh->GetDefaultOrder();

        // FSIInterface = false;
        // fSideInBoundary = -1;
        // fNeighborElements.clear();

        // IntegQuadrature nQuad(fMesh->Dimension(),fMesh->GetDefaultOrder());
        // intPointWeightFunction.resize(nQuad.getNumberOfIntegrationPoints());
        // fIntPointDistFunction.resize(nQuad.getNumberOfIntegrationPoints());
        // intPointWeightFunctionPrev.resize(nQuad.getNumberOfIntegrationPoints());

        // intPointWeightFunction.fill(1.);
        // intPointWeightFunctionPrev.fill(1.);
    
        // getIntegPointCoordinates();

    };

    ~Element() = default;

    
    virtual void ComputeElContribution(MatrixDouble &Stiffness, VecDouble &Rhs) = 0;
    virtual void ComputeElContribution(std::vector<MatrixDouble> &Stiffness, std::vector<VecDouble> &Rhs) = 0;
    virtual void ApplyBC(MatrixDouble &Stiffness, VecDouble &Rhs){};
    virtual void ApplyBC(std::vector<MatrixDouble> &Stiffness, std::vector<VecDouble> &Rhs){};
    virtual void ComputeStiffness(int &index, MatrixDouble &Stiffness){};
    virtual void ComputeStiffness(int &index, std::vector<MatrixDouble> &Stiffness){};
    virtual void ComputeResidual(int &index, VecDouble &Rhs){};
    virtual void ComputeResidual(int &index, std::vector<VecDouble> &Rhs){};
    virtual void ComputeError(VecDouble &errors){};
    

    /// Sets the element connectivity
    /// @param int* element connectivity
    void setConnectivity(VecInt &connect){fConnect = connect;};

    /// Gets the element connectivity
    /// @return element connectivity
    VecInt &getConnectivity(){return fConnect;};
    
    CompMesh* Mesh() {return fMesh;}
    
    void SetMesh(CompMesh* mesh){fMesh = mesh;}

    int64_t &Index(){return fIndex;}

    /// Compute and store the spatial jacobian matrix
    /// @param bounded_vector integration point adimensional coordinates
    virtual void ComputeJacobian(int index) = 0;
    virtual void ComputeCurrentJacobian(int index) = 0;

    /// Compute and store the shape function spatial derivatives
    /// @param bounded_vector integration point adimensional coordinates
    virtual void ComputeSpatialDerivatives() = 0;
    virtual void ComputeCurrentSpatialDerivatives() = 0;

    virtual void interpolateSolution(int &index, VecDouble &u_) = 0;
    virtual void interpolateSolution(VecDouble &phi, VecDouble &u_) = 0;
    virtual void interpolateSolDerivatives(MatrixDouble &du_dx) = 0;
    


    /// Gets the integration point global coordinates
    /// @param int integration point index @return integration point coordinates
    VecDouble getIntegPointCoordinatesValue(int index){
        VecDouble aux(2);
        aux[0] = fIntPointCoordinates(index,0); 
        aux[1] = fIntPointCoordinates(index,1); 
        return aux;
    };

    /// Sets the integration point energy weight function
    /// @param int integration point index 
    /// @param double energy weight function value
    void setIntegPointWeightFunction();
    void setIntegPointWeightFunction(int index, double val){
        intPointWeightFunction[index] = val;
    };

    /// Gets the integration point energy weight function
    /// @param int integration point index @return energy weight function value
    double &getIntegPointWeightFunction(int index) {return intPointWeightFunction[index];};
    virtual int getBoundaryGroup() = 0;
    virtual int getElement() = 0;
    virtual int getConstrain(int dir) = 0;
    virtual int getElementSide() = 0;
    virtual double getConstrainValue(int dir) = 0;
    virtual void setElement(int el) = 0;
    virtual void setElementSide(int el) = 0;
    virtual void setIntersectionParameters(VecDouble &x, VecDouble &X) = 0;
    void setFSIInterface(){FSIInterface = true;};

    /// Sets the element side in boundary
    /// @param int side in boundary
    void setElemSideInBoundary(int side){fSideInBoundary = side;};

    /// Gets the element side in boundary
    /// @return side in boundary
    int &getElemSideInBoundary(){return fSideInBoundary;};

    virtual double getJacobian() = 0;

    /// Pushs back a term of the inverse incidence, i.e., an element which
    /// contains the node
    /// @param int element
    void pushNeighborElement(int el) {
        fNeighborElements.push_back(el);
    }

    /// Gets the number of elements which contains the node
    /// @return int number of elements which contains the node
    int getNumberOfNeighborElements(){
        return fNeighborElements.size();
    }

    /// Gets an specific member of the inverse incidence
    /// @param int index @return int element of the inverse incidence
    int64_t &getNeighborElement(int i){
        return fNeighborElements[i];
    }

    void clearInverseIncidence(){
        fNeighborElements.clear();
        fNeighborElements.shrink_to_fit();
    }


    void sortEraseNeighborElements(){
        std::sort(fNeighborElements.begin(), fNeighborElements.end());
        fNeighborElements.erase(std::unique(fNeighborElements.begin(),fNeighborElements.end()), fNeighborElements.end());
    }

    virtual void setBoundaryGroup(int gr) = 0;
};

#endif