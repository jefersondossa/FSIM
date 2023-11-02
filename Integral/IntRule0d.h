#ifndef INTEGRATIONRULE0D
#define INTEGRATIONRULE0D

#include <cmath>
#include <stdio.h>
#include "IntRule.h"

/**
 @brief Integration rule associated with a point
 @ingroup integration
 */
class IntRule0d : public IntRule
{
    
public:
  
    // Default Constructor of integration rule 1D
    IntRule0d();
    
    // Constructor of integration rule 1D
    IntRule0d(int order);
    
    // Method to set polynomial order of the integration rule 1D
    virtual void SetOrder(int order) override;
    
    // Dimension of the integration rule
    virtual int Dimension() const override
    {
        return 0;
    }
    
    // Return the maximum polynomial order that can be integrated exactly
    static int gMaxOrder()
    {
        return 0;
    }

    // Return the maximum polynomial order that can be integrated exactly
    virtual int MaxOrder() const override
    {
        return gMaxOrder();
    }

};


#endif 
