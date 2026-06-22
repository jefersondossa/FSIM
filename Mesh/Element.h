#ifndef ELEMENT
#define ELEMENT

#include "DataTypes.h"
#include "IntPointData.h"
#include "WeakForm.h"
#include "GeoElement.h"
#include "Connect.h"

class GeoElement;
class CompMesh;
class Connect;
class IntPointData;

class Element{
protected:
    CompMesh                *fMesh = nullptr;
    std::vector<Connect *>  fConnect;  // element connectivity
    int                 fIndex;    //Element index

    int nLocDOF = 0;

    MatrixDouble fIntPointCoordinates;

    /// The weak form associated with the element
    WeakForm *fWeakForm = nullptr;
    IntPointData  fIntegData;

    GeoElement *fReference = nullptr;

public:
    Element(){
        fIntegData.fElement = this;
    };

    virtual ~Element() = default;

    int &NLocDOF() {
        return nLocDOF;
    }

    GeoElement *Reference() {return fReference;}

     /// Compute and store the shape function spatial derivatives
    /// @param bounded_vector integration point adimensional coordinates
    virtual void ComputeSpatialDerivatives() = 0;
    virtual void ComputeCurrentSpatialDerivatives() = 0;
    virtual void ComputeHighOrderSpatialDerivatives() = 0;

    virtual double getJacobian() = 0;
    
    virtual void ComputeElContribution(MatrixDouble &Stiffness, VecDouble &Rhs) = 0;
    virtual void ComputeElContribution(MatrixDouble &Stiffness) = 0;
    virtual void ComputeElContribution(VecDouble &Rhs) = 0;
    virtual void ComputeElContribution(std::vector<MatrixDouble> &Stiffness, std::vector<VecDouble> &Rhs) = 0;
    virtual void ComputeElContribution(std::vector<MatrixDouble> &Stiffness) = 0;
    virtual void ComputeElContribution(std::vector<VecDouble> &Rhs) = 0;
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
    virtual double IntegPointWeight(int index) = 0;
    virtual double IntegPointCoordinate(int index, int dir) = 0;

    /// Sets the element connectivity
    /// @param int* element connectivity
    virtual void setConnectivity(std::vector<Connect *> connect){fConnect = connect;};
    virtual void setConnectivity(VecInt connect);

    /// Gets the element connectivity
    /// @return element connectivity
    std::vector<Connect *> &getConnectivity(){return fConnect;};

    VecInt getConnectivityIndices(){
        VecInt indices(fConnect.size());
        for (size_t i = 0; i < fConnect.size(); i++){
            indices[i] = fConnect[i]->Index();
        }
        return indices;
    }
    
    CompMesh* Mesh() {return fMesh;}
    
    void SetMesh(CompMesh* mesh){fMesh = mesh;}

    int &Index(){return fIndex;}

    virtual void interpolateSolution(int &index, VecDouble &u_) = 0;
    virtual void interpolateSolution(VecDouble &phi, VecDouble &u_) = 0;
    virtual void interpolateSolution() = 0;
    virtual void interpolateSolDerivatives(MatrixDouble &du_dx) = 0;
    virtual void interpolateSolDerivatives(MatrixDouble &dphidx, MatrixDouble &du_dx) = 0;
    virtual void interpolateSolDerivatives() = 0;
    virtual void interpolateSolDTimeDerivatives(VecDouble &du_dt, VecDouble &du_ddt) = 0;
    virtual void interpolateSolDTimeDerivatives() = 0;

    /// Gets the integration point global coordinates
    /// @param int integration point index @return integration point coordinates
    VecDouble getIntegPointCoordinatesValue(int index){
        VecDouble aux(3);
        for (int i = 0; i < 3; i++){
            aux[i] = fIntPointCoordinates(index,i);
        }
        return aux;
    };

    MatrixDouble &getIntPointCoordinates() {return fIntPointCoordinates;};

    virtual void ComputeIntegPointCoordinates() = 0;

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

    virtual int NSides() = 0;
    virtual int NSideNodes(int iside) = 0;
    virtual int SideNodeLocIndex(int side, int node) = 0;

    double &GetIntPointDistFunction(int index){
        return  fIntegData.fDistFunction[index];
    }

    virtual IntPointData &IntegrationData() {return fIntegData;}

    virtual int Dimension() = 0;

    // Method for creating a copy of the element
    virtual Element *Clone() const = 0;

    virtual void Solution(int var, VecDouble &Sol){
        fWeakForm->Solution(fIntegData,var,Sol);
    };

    virtual void Integrate(std::vector<std::string> &varNames, std::map<std::string,VecDouble> &result){
        PanicButton();
    };

    virtual ElementType Type() = 0;

    virtual int NShapeFunctions() = 0;
};


#endif