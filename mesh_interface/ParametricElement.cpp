#include "ParametricElement.h"




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

std::vector<NodeMesh*> ParametricElement::NodeVec() const
{
    return nodes_;
}
