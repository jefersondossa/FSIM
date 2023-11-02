#ifndef INTEGRATIONRULE
#define INTEGRATIONRULE

#include <cmath>
#include <stdio.h>
#include "DataTypes.h"

/**
 @defgroup integration Integration rules
 @brief Groups classes that define integration rules
 */

/**
 @brief Root class for integration rules
 This class stores the integration points and weights an integration rule
 @ingroup integration
 */
class IntRule
{
  
protected:
    
    // Polynomial order of the integration rule
    int fOrder;
    
    // Number of integration points for this object
    MatrixDouble fPoints;
    
    // Weight of the integration point
    VecDouble fWeights;

public:
  
    // Default Constructor of integration rule
    IntRule();
    
    // Constructor of integration rule
    IntRule(int order);
    
    // Destructor of integration rule
    virtual ~IntRule();
    
    // Operator of copy
    virtual IntRule &operator=(const IntRule &copy);
    
    // Copy constructor of integration rule
    IntRule(const IntRule &copy);
    
    // Dimension of the integration rule
    virtual int Dimension() const = 0;
    
    // Return the maximum polynomial order that can be integrated exactly
    virtual int MaxOrder() const = 0;
    
    // Method to set polynomial order of the integration rule
    virtual void SetOrder(int order)
    {
        fOrder=order;
    }
    
    // Method to get polynomial order of the integration rule
    virtual int GetOrder() const
    {
        return fOrder;
    }
    
    // Method to return the number of integration points
    virtual int NPoints() const;
    
    // Function returning coordinates and weights of integration points
    virtual void Point(int p, VecDouble &co, double &weight) const;
    virtual double PointList(int i, int j) const;
    virtual double WeightList(int i) const;
    
    // Function for printing results
    virtual void Print(std::ostream &out) const;
    
};


#endif
