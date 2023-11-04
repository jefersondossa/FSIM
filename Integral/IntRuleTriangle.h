#ifndef INTEGRATIONRULETRIANGLE
#define INTEGRATIONRULETRIANGLE

#include <stdio.h>
#include "IntRule.h"

/**
@brief Integration rule associated with a triangle
@ingroup integration
*/
class IntRuleTriangle : public IntRule
{
  
    public:
  
    // Default Constructor of integration rule for triangle elements
    IntRuleTriangle();
  
    // Constructor of integration rule for triangle elements
    IntRuleTriangle(int order);
  
    // Dimension of the integration rule
    virtual int Dimension() const override{
        return 2;
    }
    
    // Return the maximum polynomial order that can be integrated exactly
    static int gMaxOrder()
    {
        return 10;
    }

    // Return the maximum polynomial order that can be integrated exactly
    virtual int MaxOrder() const override
    {
        return gMaxOrder();
    }
    // Method to set polynomial order of the integration rule for triangle elements
    virtual void SetOrder(int order) override;
  
};


#endif 
