#ifndef ELLIPSE_H
#define ELLIPSE_H

#include "Line.h"

class Ellipse : public Line
{
private:
    double rx_;
    double ry_;
    double angle_;
    
    public:
    
    Ellipse(const int& index, const std::string& name, std::vector<Point*> center, double rx, double ry, double angle = 0, const bool& discretization = true);

    ~Ellipse();

    Ellipse* operator-();

    std::string getGmshCode() override;

    Point* getMiddlePoint();
};

#endif