#ifndef INTEGRATIONRULE1D
#define INTEGRATIONRULE1D

#include <cmath>
#include <stdio.h>
#include "IntRule.h"

/**
@brief Integration rule associated with a line
@ingroup integration
*/
class IntRule1d : public IntRule
{
    
public:
  
    // Default Constructor of integration rule 1D
    IntRule1d();
    
    // Constructor of integration rule 1D
    IntRule1d(int order);
    IntRule1d(int nPoints, int order);
    
    // Method to set polynomial order of the integration rule 1D
    virtual void SetOrder(int order) override;
    
    // Dimension of the integration rule
    virtual int Dimension() const override
    {
        return 1;
    }
    
    // Return the maximum polynomial order that can be integrated exactly
    static int gMaxOrder()
    {
        return 20;
    }

    // Return the maximum polynomial order that can be integrated exactly
    virtual int MaxOrder() const override
    {
        return gMaxOrder();
    }

    // Integration rule 1D method obtained from Numerical Recipes
    static void gauleg(const REAL x1, const REAL x2, VecDouble &x, VecDouble &w);
    
};


#endif 
