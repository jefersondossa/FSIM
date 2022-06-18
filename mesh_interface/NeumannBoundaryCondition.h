#ifndef NEUMANN_BC_H
#define NEUMANN_BC_H

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


#endif

