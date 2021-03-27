#pragma once
#include "Point.h"
#include "Line.h"
#include "Surface.h"

class GeometricNeumann
{
    protected:
        int index_;
        double valueX_;
        double valueY_;
        double valueZ_;

    public:

        GeometricNeumann(const int index, const double valueX, const double valueY, const double valueZ);

        virtual ~GeometricNeumann() = 0;

        int getIndex() const;

        double getValueX() const;

        double getValueY() const;

        double getValueZ() const;

        virtual Point* getPoint() const = 0;

        virtual Line* getLine() const = 0;

        virtual Surface* getSurface() const = 0;
};

class GeometricPointNeumann : public GeometricNeumann
{
    private:
        Point* point_;

    public:

        GeometricPointNeumann(const int index, Point* const point, const double valueX, const double valueY, const double valueZ);

        ~GeometricPointNeumann() override;

        Point* getPoint() const override;

        Line* getLine() const override;

        Surface* getSurface() const override;
};

class GeometricLineNeumann : public GeometricNeumann
{
    private:
        Line* line_;

    public:

        GeometricLineNeumann(const int index, Line* const line, const double valueX, const double valueY, const double valueZ);

        ~GeometricLineNeumann() override;

        Point* getPoint() const override;

        Line* getLine() const override;

        Surface* getSurface() const override;
};

class GeometricSurfaceNeumann : public GeometricNeumann
{
    private:
        Surface* surface_;
        
    public:

        GeometricSurfaceNeumann(const int index, Surface* const surface, const double valueX, const double valueY, const double valueZ);

        ~GeometricSurfaceNeumann() override;

        Point* getPoint() const override;

        Line* getLine() const override;

        Surface* getSurface() const override;
};


GeometricNeumann::GeometricNeumann(const int index, const double valueX, const double valueY, const double valueZ)
    : index_(index), valueX_(valueX), valueY_(valueY), valueZ_(valueZ) {}

GeometricNeumann::~GeometricNeumann() {}

int GeometricNeumann::getIndex() const
{
    return index_;
}

double GeometricNeumann::getValueX() const
{
    return valueX_;
}
        
double GeometricNeumann::getValueY() const
{
    return valueY_;
}

double GeometricNeumann::getValueZ() const
{
    return valueZ_;
}

GeometricPointNeumann::GeometricPointNeumann(const int index, Point* const point, const double valueX, const double valueY, const double valueZ)
    : GeometricNeumann(index, valueX, valueY, valueZ), point_(point) {}

GeometricPointNeumann::~GeometricPointNeumann() {}

Point* GeometricPointNeumann::getPoint() const
{
    return point_;
}

Line* GeometricPointNeumann::getLine() const
{
    return nullptr;
}

Surface* GeometricPointNeumann::getSurface() const
{
    return nullptr;
}

GeometricLineNeumann::GeometricLineNeumann(const int index, Line* const line, const double valueX, const double valueY, const double valueZ)
    : GeometricNeumann(index, valueX, valueY, valueZ), line_(line) {}

GeometricLineNeumann::~GeometricLineNeumann() {}

Point* GeometricLineNeumann::getPoint() const
{
    return nullptr;
}

Line* GeometricLineNeumann::getLine() const
{
    return line_;
}

Surface* GeometricLineNeumann::getSurface() const
{
    return nullptr;
}

GeometricSurfaceNeumann::GeometricSurfaceNeumann(const int index, Surface* const surface, const double valueX, const double valueY, const double valueZ)
    : GeometricNeumann(index, valueX, valueY, valueZ), surface_(surface) {}

GeometricSurfaceNeumann::~GeometricSurfaceNeumann() {}

Point* GeometricSurfaceNeumann::getPoint() const
{
    return nullptr;
}

Line* GeometricSurfaceNeumann::getLine() const
{
    return nullptr;
}

Surface* GeometricSurfaceNeumann::getSurface() const
{
    return surface_;
}
