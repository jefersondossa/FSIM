#ifndef INTEGRATIONRULEHEXA
#define INTEGRATIONRULEHEXA

#include <stdio.h>
#include "IntRule.h"

/**
@brief Integration rule associated with a tetrahedron
@ingroup integration
*/
class IntRuleHexahedron : public IntRule
{
  
    public:
    
    // Default Constructor of integration rule for tetrahedron elements
    IntRuleHexahedron();
    
    // Constructor of integration rule for tetrahedron elements
    IntRuleHexahedron(int order);

    // Dimension of the integration rule
    virtual int Dimension() const override{
        return 3;
    }
    
    // Return the maximum polynomial order that can be integrated exactly
    static int gMaxOrder()
    {
        return 14;
    }

    // Return the maximum polynomial order that can be integrated exactly
    virtual int MaxOrder() const override
    {
        return gMaxOrder();
    }
    // Method to set polynomial order of the integration rule for tetrahedro elements
    virtual void SetOrder(int order) override;
    
};


#endif 
