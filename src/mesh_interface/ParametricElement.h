#pragma once
#include "NodeMesh.h"
// #include "Quadratures.h"
#include <vector>
#include <cmath>

enum PartitionOfUnity
{
    //point
    VERTEX,
    // line partition of unity
    L2, L3, L4,
    // triangular partition of unity
    T3, T6, T10,
    // quadrilateral partition of unity
    Q4, Q9, Q16,
    // tetrahedron partition of unity
    TET4, TET10, TET20
};

class ParametricElement
{
    public:
        ParametricElement(const int& index, const std::vector<NodeMesh*>& nodes);

        virtual ~ParametricElement();

        void setIndex(const int& index);

        void setOrder(const int& order);

        void setNumberOfIntegrationPoints(const int& numberOfIntegrationPoints);

        void setPlot(const bool& plot);

        void setPartitionOfUnity(const PartitionOfUnity& partitionOfUnity);

        void setNodes(const std::vector<NodeMesh*>& nodes);

        int getIndex() const;

        int getOrder() const;

        int getNumberOfIntegrationPoints() const;

        bool getPlot() const;

        int getNumberOfNodes() const;

        PartitionOfUnity getPartitionOfUnity() const;

        std::vector<NodeMesh*> getNodes() const;

        virtual std::vector<NodeMesh*> getBoundaryNodes(const int& boundaryIndex) const = 0;
    
    protected:
        int index_;
        int order_;
        int numberOfIntegrationPoints_;
        bool plot_;
        std::vector<NodeMesh*> nodes_;
        PartitionOfUnity partitionOfUnity_;
};















ParametricElement::ParametricElement(const int& index, const std::vector<NodeMesh*>& nodes)
    : index_(index), nodes_(nodes), plot_(false) {}

ParametricElement::~ParametricElement() {}

void ParametricElement::setIndex(const int& index)
{
    index_ = index;
}

void ParametricElement::setOrder(const int& order)
{
    order_ = order;
}

void ParametricElement::setNumberOfIntegrationPoints(const int& numberOfIntegrationPoints)
{
    numberOfIntegrationPoints_ = numberOfIntegrationPoints;
}

void ParametricElement::setPlot(const bool& plot)
{
    plot_ = plot;
}

void ParametricElement::setPartitionOfUnity(const PartitionOfUnity& partitionOfUnity)
{
    partitionOfUnity_ = partitionOfUnity;
}

void ParametricElement::setNodes(const std::vector<NodeMesh *>& nodes)
{
    nodes_ = nodes;
}

int ParametricElement::getIndex() const
{
    return index_;
}

int ParametricElement::getOrder() const
{
    return order_;
}

int ParametricElement::getNumberOfIntegrationPoints() const
{
    return numberOfIntegrationPoints_;
}

bool ParametricElement::getPlot() const
{
    return plot_;
}

int ParametricElement::getNumberOfNodes() const
{
    return nodes_.size();
}

PartitionOfUnity ParametricElement::getPartitionOfUnity() const
{
    return partitionOfUnity_;
}

std::vector<NodeMesh*> ParametricElement::getNodes() const
{
    return nodes_;
}
