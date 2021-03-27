#pragma once
#include "NodeMesh.h"
// #include "ParametricLineElement.h"
// #include "ParametricSurfaceElement.h"
#include <cmath>
#include <iostream>

class NeumannBoundaryCondition
{
    protected:
        int index_;
        int ndofs_;
        double forces_[3];

    public:

        NeumannBoundaryCondition(const int index, const int ndofs, const double valueX, const double valueY, const double valueZ);

        virtual ~NeumannBoundaryCondition() = 0;

        int getNumberOfDOFs() const;

        double getForce(const int& dof) const;

        virtual void getNodalForce(const int& dimension, int*& indexes, double*& values) const = 0;
};

class PointLoad : public NeumannBoundaryCondition
{
    private:
        NodeMesh* node_;
    
    public:

        PointLoad(const int index, const int ndofs, NodeMesh* const node, const double valueX, const double valueY, const double valueZ);

        ~PointLoad() override;

        void getNodalForce(const int& dimension, int*& indexes, double*& values) const override;
};

// class LineLoad : public NeumannBoundaryCondition
// {
//     private:
//         ParametricLineElement* element_;

//     public:

//         LineLoad(const int index, const int ndofs, ParametricLineElement* const element, const double valueX, const double valueY, const double valueZ);

//         ~LineLoad() override;

//         void getNodalForce(const int& dimension, int*& indexes, double*& values) const override;
// };

// class SurfaceLoad : public NeumannBoundaryCondition
// {
//     private:
//         ParametricSurfaceElement* element_;

//     public:

//         SurfaceLoad(const int index, const int ndofs, ParametricSurfaceElement* const element, const double valueX, const double valueY, const double valueZ);

//         ~SurfaceLoad() override;
        
//         void getNodalForce(const int& dimension, int*& indexes, double*& values) const override;
// };








NeumannBoundaryCondition::NeumannBoundaryCondition(const int index, const int ndofs, const double valueX, const double valueY, const double valueZ)
    : index_(index), ndofs_(ndofs), forces_{valueX, valueY, valueZ}  {}

NeumannBoundaryCondition::~NeumannBoundaryCondition() {}

int NeumannBoundaryCondition::getNumberOfDOFs() const
{
    return ndofs_;
}

double NeumannBoundaryCondition::getForce(const int& dof) const
{
    return forces_[dof];
}

// PointLoad::PointLoad(const int index, const int ndofs, NodeMesh* const node, const double valueX, const double valueY, const double valueZ)
//     : NeumannBoundaryCondition(index, ndofs, valueX, valueY, valueZ), node_(node) {}

// PointLoad::~PointLoad() {}

// void PointLoad::getNodalForce(const int& dimension, int*& indexes, double*& values) const
// {
//     indexes = new int[ndofs_];
//     values = new double[ndofs_];
//     for (int i = 0; i < dimension; i++)
//     {
//         indexes[i] = node_->getDegreeOfFreedom(i)->getIndex();
//         values[i] = forces_[i];
//     }
// }

// LineLoad::LineLoad(const int index, const int ndofs, ParametricLineElement* const element, const double valueX, const double valueY, const double valueZ)
//     : NeumannBoundaryCondition(index, ndofs, valueX, valueY, valueZ), element_(element) {}

// LineLoad::~LineLoad() {}

// void LineLoad::getNodalForce(const int& dimension, int*& indexes, double*& values) const
// {
//     indexes = new int[ndofs_];
//     values = new double[ndofs_];

//     for (int i = 0; i < ndofs_; ++i)
//         values[i] = 0.0;
    
//     const int numberOfNodes = element_->getNumberOfNodes();
//     const std::vector<NodeMesh*>& nodes = element_->getNodes();
    
//     const int numberOfIntegrationPoints = element_->getNumberOfIntegrationPoints();
//     double *xsi, *weight;
//     quadratures::lineQuadrature(numberOfIntegrationPoints, xsi, weight);

//     // integration loop
//     for (int ip = 0; ip < numberOfIntegrationPoints; ip++)
//     {
//         double *phi, *dphi_dxsi;
//         element_->getShapeFunctions(xsi[ip], phi, dphi_dxsi);
        
//         double dx_dxsi[dimension];
//         for (int i = 0; i < dimension; i++)
//             dx_dxsi[i] = 0.0;
        
//         double interpolatedForce[dimension];
//         for (int i = 0; i < dimension; i++)
//             interpolatedForce[i] = 0.0;  
    
//         for (int i = 0; i < numberOfNodes; ++i)
//         {
//             double coord;
//             for (int j = 0; j < dimension; j++)
//             {
//                 coord = nodes[i]->getDegreeOfFreedom(j)->getInitialValue();
//                 dx_dxsi[j] += dphi_dxsi[i] * coord;
//                 interpolatedForce[j] += phi[i] * forces_[j];
//             }
//         }
//         double j0 = 0;
//         for (int i = 0; i < dimension; i++)
//         {
//             j0 += dx_dxsi[i] * dx_dxsi[i];
//         }
//         j0 = sqrt(j0);

//         for (int i = 0; i < numberOfNodes; i++)
//         {
//             for (int j = 0; j < dimension; j++)
//             {
//                 values[dimension*i+j] += interpolatedForce[j] * phi[i] * weight[ip] * j0;
//             }
//         }
//         delete[] phi;
//         delete[] dphi_dxsi;
//     }
//     delete[] xsi;
//     delete[] weight;

//     for (int i = 0; i < numberOfNodes; i++)
//         for (int j = 0; j < dimension; j++)
//             indexes[dimension*i+j] = nodes[i]->getDegreeOfFreedom(j)->getIndex();
// }

// SurfaceLoad::SurfaceLoad(const int index, const int ndofs, ParametricSurfaceElement* const element, const double valueX, const double valueY, const double valueZ)
//     : NeumannBoundaryCondition(index, ndofs, valueX, valueY, valueZ), element_(element) {}

// SurfaceLoad::~SurfaceLoad() {}

// void SurfaceLoad::getNodalForce(const int& dimension, int*& indexes, double*& values) const
// {
//     indexes = new int[ndofs_];
//     values = new double[ndofs_];

//     for (int i = 0; i < ndofs_; i++)
//         values[i] = 0.0;
    
//     const int numberOfNodes = element_->getNumberOfNodes();
//     const std::vector<NodeMesh*> nodes = element_->getNodes();
    
//     const int numberOfIntegrationPoints = element_->getNumberOfIntegrationPoints();
//     double *xsi1, *xsi2, *weight;
//     quadratures::triangleQuadrature(numberOfIntegrationPoints, xsi1, xsi2, weight);
//     // integration loop
//     for (int ip = 0; ip < numberOfIntegrationPoints; ip++)
//     {
//         double *phi, *dphi_dxsi1, *dphi_dxsi2;
//         element_->getShapeFunctions(xsi1[ip], xsi2[ip], phi, dphi_dxsi1, dphi_dxsi2);
//         double dx_dxsi1[3], dx_dxsi2[3];
//         dx_dxsi1[0] = 0.0; dx_dxsi1[1] = 0.0; dx_dxsi1[2] = 0.0;
//         dx_dxsi2[0] = 0.0; dx_dxsi2[1] = 0.0; dx_dxsi2[2] = 0.0;
//         double interpolatedForce[3];
//         interpolatedForce[0] = 0.0; interpolatedForce[1] = 0.0; interpolatedForce[2] = 0.0;
        
//         for (int i = 0; i < numberOfNodes; i++)
//         {
//             double coords[3];
//             coords[0] = nodes[i]->getDegreeOfFreedom(0)->getInitialValue();
//             coords[1] = nodes[i]->getDegreeOfFreedom(1)->getInitialValue();
//             coords[2] = nodes[i]->getDegreeOfFreedom(2)->getInitialValue();

//             dx_dxsi1[0] += dphi_dxsi1[i] * coords[0];
//             dx_dxsi1[1] += dphi_dxsi1[i] * coords[1];
//             dx_dxsi1[2] += dphi_dxsi1[i] * coords[2];
            
//             dx_dxsi2[0] += dphi_dxsi2[i] * coords[0];
//             dx_dxsi2[1] += dphi_dxsi2[i] * coords[1];
//             dx_dxsi2[2] += dphi_dxsi2[i] * coords[2];

//             interpolatedForce[0] += phi[i] * forces_[0];
//             interpolatedForce[1] += phi[i] * forces_[1];
//             interpolatedForce[2] += phi[i] * forces_[2];
//         }
//         double j0 = (dx_dxsi1[0]*dx_dxsi1[0] + dx_dxsi1[1]*dx_dxsi1[1] + dx_dxsi1[2]*dx_dxsi1[2]) * 
//                     (dx_dxsi2[0]*dx_dxsi2[0] + dx_dxsi2[1]*dx_dxsi2[1] + dx_dxsi2[2]*dx_dxsi2[2]) - 
// 					(dx_dxsi1[0]*dx_dxsi2[0] + dx_dxsi1[1]*dx_dxsi2[1] + dx_dxsi1[2]*dx_dxsi2[2]) * 
//                     (dx_dxsi2[0]*dx_dxsi1[0] + dx_dxsi2[1]*dx_dxsi1[1] + dx_dxsi2[2]*dx_dxsi1[2]);
// 		j0 = sqrt(j0);

//         for (int i = 0; i < numberOfNodes; i++)
//         {
//             values[3*i+0] += interpolatedForce[0] * phi[i] * weight[ip] * j0;
//             values[3*i+1] += interpolatedForce[1] * phi[i] * weight[ip] * j0;
//             values[3*i+2] += interpolatedForce[2] * phi[i] * weight[ip] * j0;
//         }
//         delete[] phi;
//         delete[] dphi_dxsi1;
//         delete[] dphi_dxsi2;
//     }
//     delete[] xsi1;
//     delete[] xsi2;
//     delete[] weight;

//     for (int i = 0; i < numberOfNodes; i++)
//     {
//         indexes[3*i+0] = nodes[i]->getDegreeOfFreedom(0)->getIndex();
//         indexes[3*i+1] = nodes[i]->getDegreeOfFreedom(1)->getIndex();
//         indexes[3*i+2] = nodes[i]->getDegreeOfFreedom(2)->getIndex();
//     }
// }





