#ifndef ELEMENT
#define ELEMENT

#include "DataTypes.h"
#include "CompMesh.h"
#include "IntPointData.h"
#include "WeakForm.h"

class CompMesh;

class Element{
protected:
    CompMesh *fMesh;
    VecInt        fConnect; //Velocity mesh connectivity 
    int64_t       fIndex;             //Element index

    int nLocDOF = 0;

    MatrixDouble fIntPointCoordinates;
    int DEG;

    VecDouble     xK, XK;

    std::vector<int64_t> fNeighborElements;
    /// The weak form associated with the element
    WeakForm *fWeakForm;
    IntPointData  fIntegData;

    int fPrintType;
public:

    bool          FSIInterface;  

public:
    Element() = default;

    // Element(int64_t index, VecInt &connect, CompMesh* mesh){
        
    //     fMesh = mesh;
    //     fConnect.resize(connect.size());
    //     fIndex = index;
    //     for (int i = connect.size(); i--; ) fConnect[i] = connect[i];

    // };

    ~Element() = default;
    int &NLocDOF() {return nLocDOF;}
    virtual const int &NElNodes() = 0;

    int &PrintType() {return fPrintType;}
    
    virtual void ComputeElContribution(MatrixDouble &Stiffness, VecDouble &Rhs) = 0;
    virtual void ComputeElContribution(std::vector<MatrixDouble> &Stiffness, std::vector<VecDouble> &Rhs) = 0;
    virtual void ApplyBC(MatrixDouble &Stiffness, VecDouble &Rhs){};
    virtual void ApplyBC(std::vector<MatrixDouble> &Stiffness, std::vector<VecDouble> &Rhs){};
    virtual void ComputeStiffness(int &index, MatrixDouble &Stiffness){};
    virtual void ComputeStiffness(int &index, std::vector<MatrixDouble> &Stiffness){};
    virtual void ComputeResidual(int &index, VecDouble &Rhs){};
    virtual void ComputeResidual(int &index, std::vector<VecDouble> &Rhs){};
    virtual void ComputeError(VecDouble &errors){};
    
    virtual double InterpolateVariable(VecDouble &nValues, int point) = 0;

    WeakForm *GetWeakForm(){return fWeakForm;}
    void SetWeakForm(WeakForm *weak){fWeakForm=weak;}

    virtual void ComputeIntPointDistFunction(VecDouble &nodalval) = 0;
    virtual int getNumberOfIntegrationPoints() = 0;

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
    virtual void ComputeJacobian() = 0;
    virtual void ComputeCurrentJacobian() = 0;

    /// Compute and store the shape function spatial derivatives
    /// @param bounded_vector integration point adimensional coordinates
    virtual void ComputeSpatialDerivatives() = 0;
    virtual void ComputeCurrentSpatialDerivatives() = 0;

    virtual void interpolateSolution(int &index, VecDouble &u_) = 0;
    virtual void interpolateSolution(VecDouble &phi, VecDouble &u_) = 0;
    virtual void interpolateSolution() = 0;
    virtual void interpolateSolDerivatives(MatrixDouble &du_dx) = 0;
    virtual void interpolateSolDerivatives(MatrixDouble &dphidx, MatrixDouble &du_dx) = 0;
    virtual void interpolateSolDerivatives() = 0;
    


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
    virtual void setIntegPointWeightFunction() = 0;
    void setIntegPointWeightFunction(int index, double val){
        fIntegData.fWeightFunction[index] = val;
    };

    /// Gets the integration point energy weight function
    /// @param int integration point index @return energy weight function value
    double &getIntegPointWeightFunction(int index) {return fIntegData.fWeightFunction[index];};
    virtual int getElement() = 0;
    virtual int getElementSide() = 0;
    virtual void setElement(int el) = 0;
    virtual void setElementSide(int el) = 0;
    virtual void setIntersectionParameters(VecDouble &x, VecDouble &X) = 0;
    void setFSIInterface(){FSIInterface = true;};

    /// Gets the coordinates intersection parameters
    /// @return minimum and maximum coordinates
    std::pair<VecDouble,VecDouble> getXIntersectionParameter() {return std::make_pair(xK,XK);};

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

    virtual int NCornerNodes() = 0;

    double &GetIntPointDistFunction(int index){
        return  fIntegData.fDistFunction[index];
    }

    IntPointData &IntegrationData() {return fIntegData;}

    virtual int Dimension() = 0;

    // Method for creating a copy of the element
    virtual Element *Clone() const = 0;

    virtual VecDouble NodeCoord(int inode) = 0;

    void Solution(int var, VecDouble &Sol){
        fWeakForm->Solution(fIntegData,var,Sol);
    };
};


#endif