#ifndef PAR_ELEMENT_H
#define PAR_ELEMENT_H

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

#endif











