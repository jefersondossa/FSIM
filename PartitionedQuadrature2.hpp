//------------------------------------------------------------------------------
// 
//                   Jeferson W D Fernandes and Rodolfo A K Sanches
//                             University of Sao Paulo
//                           (C) 2017 All Rights Reserved
//
// <LicenseText>
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//------------------------------QUADRATURE POINTS-------------------------------
//------------------------------------------------------------------------------

#ifndef PART_QUADRATURE_H
#define PART_QUADRATURE_H

#include "QuadraticShapeFunction.hpp"
#include "LinearShapeFunction.hpp"

template<int DIM>
class PartQuadrature{
public:
    //Defines the type "PointCoord" which stores the 
    //integration points coordinates
    typedef ublas::bounded_matrix<double, 112,DIM>  PointCoord;

    //Defines the type "PointWeight" which stores the
    //integration points weights
    typedef ublas::bounded_vector<double, 112>      PointWeight;

    //Integration point logical vector
    typedef ublas::bounded_vector<bool, 112>        PointLogical;

    //Defines the numerical integration iterator
    typedef typename PointWeight::iterator              QuadratureListIt;

    //Defines vector of nodal values
    typedef ublas::bounded_vector<double, 4*DIM-2>      NodalValuesQuad;
    typedef ublas::bounded_vector<double, DIM+1>        NodalValuesLin;

public:
    //Returns the index of the first integration point
    QuadratureListIt begin() {
        return pointWeight.begin();
    }

    //Returns the index of the last integration point
    QuadratureListIt end() {
        return pointWeight.end();
    }

    //Returns the integration point coordinate
    double PointList(int i, int j); 

    //Retuns the integration point weight
    double WeightList(int i);
  
    //Interpolate variables
    double interpolateQuadraticVariable(NodalValuesQuad nValues, int point);
    double interpolateLinearVariable(NodalValuesLin nValues, int point);

private:
    //List of integration points coordinates
    PointCoord pointCoord;

    //List of integration points weights
    PointWeight pointWeight;

    //Defines shape functions
    QuadShapeFunction<DIM> shapeQuad;
    LinShapeFunction<DIM>  shapeLin;

    //Values of velocity shape functins
    typename QuadShapeFunction<DIM>::Values      phi_;     
    //Values of pressure shape functins
    typename LinShapeFunction<DIM>::Values       phip_;    
};

//------------------------------------------------------------------------------
//--------------------------------IMPLEMENTATION--------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//-----------------------QUADRATURE POINTS - COORDINATES------------------------
//------------------------------------------------------------------------------
template<>
double PartQuadrature<2>::PointList(int i, int j){

pointCoord(	0	,0) = 	0.083333333333333 	;
pointCoord(	1	,0) = 	0.025321626830864 	;
pointCoord(	2	,0) = 	0.199356746338272 	;
pointCoord(	3	,0) = 	0.025321626830864 	;
pointCoord(	4	,0) = 	0.117535516026279 	;
pointCoord(	5	,0) = 	0.117535516026279 	;
pointCoord(	6	,0) = 	0.014928967947443 	;
pointCoord(	7	,0) = 	0.166666666666667 	;
pointCoord(	8	,0) = 	0.050643253661728 	;
pointCoord(	9	,0) = 	0.224678373169136 	;
pointCoord(	10	,0) = 	0.224678373169136 	;
pointCoord(	11	,0) = 	0.132464483973721 	;
pointCoord(	12	,0) = 	0.235071032052557 	;
pointCoord(	13	,0) = 	0.132464483973721 	;
pointCoord(	14	,0) = 	0.333333333333333 	;
pointCoord(	15	,0) = 	0.275321626830864 	;
pointCoord(	16	,0) = 	0.449356746338272 	;
pointCoord(	17	,0) = 	0.275321626830864 	;
pointCoord(	18	,0) = 	0.367535516026279 	;
pointCoord(	19	,0) = 	0.367535516026279 	;
pointCoord(	20	,0) = 	0.264928967947442 	;
pointCoord(	21	,0) = 	0.416666666666667 	;
pointCoord(	22	,0) = 	0.300643253661728 	;
pointCoord(	23	,0) = 	0.474678373169136 	;
pointCoord(	24	,0) = 	0.474678373169136 	;
pointCoord(	25	,0) = 	0.382464483973721 	;
pointCoord(	26	,0) = 	0.485071032052558 	;
pointCoord(	27	,0) = 	0.382464483973721 	;
pointCoord(	28	,0) = 	0.583333333333333 	;
pointCoord(	29	,0) = 	0.525321626830864 	;
pointCoord(	30	,0) = 	0.699356746338272 	;
pointCoord(	31	,0) = 	0.525321626830864 	;
pointCoord(	32	,0) = 	0.617535516026279 	;
pointCoord(	33	,0) = 	0.617535516026279 	;
pointCoord(	34	,0) = 	0.514928967947442 	;
pointCoord(	35	,0) = 	0.666666666666667 	;
pointCoord(	36	,0) = 	0.550643253661728 	;
pointCoord(	37	,0) = 	0.724678373169136 	;
pointCoord(	38	,0) = 	0.724678373169136 	;
pointCoord(	39	,0) = 	0.632464483973721 	;
pointCoord(	40	,0) = 	0.735071032052558 	;
pointCoord(	41	,0) = 	0.632464483973721 	;
pointCoord(	42	,0) = 	0.833333333333333 	;
pointCoord(	43	,0) = 	0.775321626830864 	;
pointCoord(	44	,0) = 	0.949356746338272 	;
pointCoord(	45	,0) = 	0.775321626830864 	;
pointCoord(	46	,0) = 	0.867535516026279 	;
pointCoord(	47	,0) = 	0.867535516026279 	;
pointCoord(	48	,0) = 	0.764928967947442 	;
pointCoord(	49	,0) = 	0.083333333333333 	;
pointCoord(	50	,0) = 	0.025321626830864 	;
pointCoord(	51	,0) = 	0.199356746338272 	;
pointCoord(	52	,0) = 	0.025321626830864 	;
pointCoord(	53	,0) = 	0.117535516026279 	;
pointCoord(	54	,0) = 	0.117535516026279 	;
pointCoord(	55	,0) = 	0.014928967947443 	;
pointCoord(	56	,0) = 	0.166666666666667 	;
pointCoord(	57	,0) = 	0.050643253661728 	;
pointCoord(	58	,0) = 	0.224678373169136 	;
pointCoord(	59	,0) = 	0.224678373169136 	;
pointCoord(	60	,0) = 	0.132464483973721 	;
pointCoord(	61	,0) = 	0.235071032052557 	;
pointCoord(	62	,0) = 	0.132464483973721 	;
pointCoord(	63	,0) = 	0.333333333333333 	;
pointCoord(	64	,0) = 	0.275321626830864 	;
pointCoord(	65	,0) = 	0.449356746338272 	;
pointCoord(	66	,0) = 	0.275321626830864 	;
pointCoord(	67	,0) = 	0.367535516026279 	;
pointCoord(	68	,0) = 	0.367535516026279 	;
pointCoord(	69	,0) = 	0.264928967947442 	;
pointCoord(	70	,0) = 	0.416666666666667 	;
pointCoord(	71	,0) = 	0.300643253661728 	;
pointCoord(	72	,0) = 	0.474678373169136 	;
pointCoord(	73	,0) = 	0.474678373169136 	;
pointCoord(	74	,0) = 	0.382464483973721 	;
pointCoord(	75	,0) = 	0.485071032052558 	;
pointCoord(	76	,0) = 	0.382464483973721 	;
pointCoord(	77	,0) = 	0.583333333333333 	;
pointCoord(	78	,0) = 	0.525321626830864 	;
pointCoord(	79	,0) = 	0.699356746338272 	;
pointCoord(	80	,0) = 	0.525321626830864 	;
pointCoord(	81	,0) = 	0.617535516026279 	;
pointCoord(	82	,0) = 	0.617535516026279 	;
pointCoord(	83	,0) = 	0.514928967947442 	;
pointCoord(	84	,0) = 	0.083333333333333 	;
pointCoord(	85	,0) = 	0.025321626830864 	;
pointCoord(	86	,0) = 	0.199356746338272 	;
pointCoord(	87	,0) = 	0.025321626830864 	;
pointCoord(	88	,0) = 	0.117535516026279 	;
pointCoord(	89	,0) = 	0.117535516026279 	;
pointCoord(	90	,0) = 	0.014928967947443 	;
pointCoord(	91	,0) = 	0.166666666666667 	;
pointCoord(	92	,0) = 	0.050643253661728 	;
pointCoord(	93	,0) = 	0.224678373169136 	;
pointCoord(	94	,0) = 	0.224678373169136 	;
pointCoord(	95	,0) = 	0.132464483973721 	;
pointCoord(	96	,0) = 	0.235071032052557 	;
pointCoord(	97	,0) = 	0.132464483973721 	;
pointCoord(	98	,0) = 	0.333333333333333 	;
pointCoord(	99	,0) = 	0.275321626830864 	;
pointCoord(	100	,0) = 	0.449356746338272 	;
pointCoord(	101	,0) = 	0.275321626830864 	;
pointCoord(	102	,0) = 	0.367535516026279 	;
pointCoord(	103	,0) = 	0.367535516026279 	;
pointCoord(	104	,0) = 	0.264928967947442 	;
pointCoord(	105	,0) = 	0.083333333333333 	;
pointCoord(	106	,0) = 	0.025321626830864 	;
pointCoord(	107	,0) = 	0.199356746338272 	;
pointCoord(	108	,0) = 	0.025321626830864 	;
pointCoord(	109	,0) = 	0.117535516026279 	;
pointCoord(	110	,0) = 	0.117535516026279 	;
pointCoord(	111	,0) = 	0.014928967947443 	;
  
pointCoord(	0	,1) = 	0.166666666666667 	;
pointCoord(	1	,1) = 	0.050643253661728 	;
pointCoord(	2	,1) = 	0.224678373169136 	;
pointCoord(	3	,1) = 	0.224678373169136 	;
pointCoord(	4	,1) = 	0.132464483973721 	;
pointCoord(	5	,1) = 	0.235071032052557 	;
pointCoord(	6	,1) = 	0.132464483973721 	;
pointCoord(	7	,1) = 	0.083333333333333 	;
pointCoord(	8	,1) = 	0.025321626830864 	;
pointCoord(	9	,1) = 	0.025321626830864 	;
pointCoord(	10	,1) = 	0.199356746338272 	;
pointCoord(	11	,1) = 	0.014928967947443 	;
pointCoord(	12	,1) = 	0.117535516026279 	;
pointCoord(	13	,1) = 	0.117535516026279 	;
pointCoord(	14	,1) = 	0.083333333333334 	;
pointCoord(	15	,1) = 	0.025321626830864 	;
pointCoord(	16	,1) = 	0.025321626830864 	;
pointCoord(	17	,1) = 	0.199356746338272 	;
pointCoord(	18	,1) = 	0.014928967947443 	;
pointCoord(	19	,1) = 	0.117535516026279 	;
pointCoord(	20	,1) = 	0.117535516026279 	;
pointCoord(	21	,1) = 	0.166666666666667 	;
pointCoord(	22	,1) = 	0.224678373169136 	;
pointCoord(	23	,1) = 	0.050643253661728 	;
pointCoord(	24	,1) = 	0.224678373169136 	;
pointCoord(	25	,1) = 	0.132464483973721 	;
pointCoord(	26	,1) = 	0.132464483973721 	;
pointCoord(	27	,1) = 	0.235071032052557 	;
pointCoord(	28	,1) = 	0.166666666666667 	;
pointCoord(	29	,1) = 	0.050643253661728 	;
pointCoord(	30	,1) = 	0.224678373169136 	;
pointCoord(	31	,1) = 	0.224678373169136 	;
pointCoord(	32	,1) = 	0.132464483973721 	;
pointCoord(	33	,1) = 	0.235071032052557 	;
pointCoord(	34	,1) = 	0.132464483973721 	;
pointCoord(	35	,1) = 	0.083333333333334 	;
pointCoord(	36	,1) = 	0.025321626830864 	;
pointCoord(	37	,1) = 	0.025321626830864 	;
pointCoord(	38	,1) = 	0.199356746338272 	;
pointCoord(	39	,1) = 	0.014928967947443 	;
pointCoord(	40	,1) = 	0.117535516026279 	;
pointCoord(	41	,1) = 	0.117535516026279 	;
pointCoord(	42	,1) = 	0.083333333333334 	;
pointCoord(	43	,1) = 	0.025321626830864 	;
pointCoord(	44	,1) = 	0.025321626830864 	;
pointCoord(	45	,1) = 	0.199356746338272 	;
pointCoord(	46	,1) = 	0.014928967947443 	;
pointCoord(	47	,1) = 	0.117535516026279 	;
pointCoord(	48	,1) = 	0.117535516026279 	;
pointCoord(	49	,1) = 	0.333333333333333 	;
pointCoord(	50	,1) = 	0.275321626830864 	;
pointCoord(	51	,1) = 	0.275321626830864 	;
pointCoord(	52	,1) = 	0.449356746338272 	;
pointCoord(	53	,1) = 	0.264928967947442 	;
pointCoord(	54	,1) = 	0.367535516026279 	;
pointCoord(	55	,1) = 	0.367535516026279 	;
pointCoord(	56	,1) = 	0.416666666666667 	;
pointCoord(	57	,1) = 	0.474678373169136 	;
pointCoord(	58	,1) = 	0.300643253661728 	;
pointCoord(	59	,1) = 	0.474678373169136 	;
pointCoord(	60	,1) = 	0.382464483973721 	;
pointCoord(	61	,1) = 	0.382464483973721 	;
pointCoord(	62	,1) = 	0.485071032052557 	;
pointCoord(	63	,1) = 	0.416666666666667 	;
pointCoord(	64	,1) = 	0.300643253661728 	;
pointCoord(	65	,1) = 	0.474678373169136 	;
pointCoord(	66	,1) = 	0.474678373169136 	;
pointCoord(	67	,1) = 	0.382464483973721 	;
pointCoord(	68	,1) = 	0.485071032052558 	;
pointCoord(	69	,1) = 	0.382464483973721 	;
pointCoord(	70	,1) = 	0.333333333333333 	;
pointCoord(	71	,1) = 	0.275321626830864 	;
pointCoord(	72	,1) = 	0.275321626830864 	;
pointCoord(	73	,1) = 	0.449356746338272 	;
pointCoord(	74	,1) = 	0.264928967947442 	;
pointCoord(	75	,1) = 	0.367535516026279 	;
pointCoord(	76	,1) = 	0.367535516026279 	;
pointCoord(	77	,1) = 	0.333333333333333 	;
pointCoord(	78	,1) = 	0.275321626830864 	;
pointCoord(	79	,1) = 	0.275321626830864 	;
pointCoord(	80	,1) = 	0.449356746338272 	;
pointCoord(	81	,1) = 	0.264928967947442 	;
pointCoord(	82	,1) = 	0.367535516026279 	;
pointCoord(	83	,1) = 	0.367535516026279 	;
pointCoord(	84	,1) = 	0.666666666666667 	;
pointCoord(	85	,1) = 	0.550643253661728 	;
pointCoord(	86	,1) = 	0.724678373169136 	;
pointCoord(	87	,1) = 	0.724678373169136 	;
pointCoord(	88	,1) = 	0.632464483973721 	;
pointCoord(	89	,1) = 	0.735071032052558 	;
pointCoord(	90	,1) = 	0.632464483973721 	;
pointCoord(	91	,1) = 	0.583333333333333 	;
pointCoord(	92	,1) = 	0.525321626830864 	;
pointCoord(	93	,1) = 	0.525321626830864 	;
pointCoord(	94	,1) = 	0.699356746338272 	;
pointCoord(	95	,1) = 	0.514928967947442 	;
pointCoord(	96	,1) = 	0.617535516026279 	;
pointCoord(	97	,1) = 	0.617535516026279 	;
pointCoord(	98	,1) = 	0.583333333333333 	;
pointCoord(	99	,1) = 	0.525321626830864 	;
pointCoord(	100	,1) = 	0.525321626830864 	;
pointCoord(	101	,1) = 	0.699356746338272 	;
pointCoord(	102	,1) = 	0.514928967947442 	;
pointCoord(	103	,1) = 	0.617535516026279 	;
pointCoord(	104	,1) = 	0.617535516026279 	;
pointCoord(	105	,1) = 	0.833333333333333 	;
pointCoord(	106	,1) = 	0.775321626830864 	;
pointCoord(	107	,1) = 	0.775321626830864 	;
pointCoord(	108	,1) = 	0.949356746338272 	;
pointCoord(	109	,1) = 	0.764928967947442 	;
pointCoord(	110	,1) = 	0.867535516026279 	;
pointCoord(	111	,1) = 	0.867535516026279 	;

 
    return pointCoord(i,j);
};

template<>
double PartQuadrature<3>::PointList(int i, int j){
    
    const double a = (1. + sqrt(5. / 14.)) / 4.;
    const double b = (1. - sqrt(5. / 14.)) / 4.;

    pointCoord(0,0) = 1. / 4.;
    pointCoord(0,1) = 1. / 4.;
    pointCoord(0,2) = 1. / 4.;

    pointCoord(1,0) = 11. / 14.;
    pointCoord(1,1) = 1. / 14.;
    pointCoord(1,2) = 1. / 14.;

    pointCoord(2,0) = 1. / 14.;
    pointCoord(2,1) = 11. / 14.;
    pointCoord(2,2) = 1. / 14.;

    pointCoord(3,0) = 1. / 14.;
    pointCoord(3,1) = 1. / 14.;
    pointCoord(3,2) = 11. / 14.;

    pointCoord(4,0) = 1. / 14.;
    pointCoord(4,1) = 1. / 14.;
    pointCoord(4,2) = 1. / 14.;

    pointCoord(5,0) = a;
    pointCoord(5,1) = a;
    pointCoord(5,2) = b;

    pointCoord(6,0) = a;
    pointCoord(6,1) = b;
    pointCoord(6,2) = a;

    pointCoord(7,0) = a;
    pointCoord(7,1) = b;
    pointCoord(7,2) = b;

    pointCoord(8,0) = b;
    pointCoord(8,1) = a;
    pointCoord(8,2) = a;

    pointCoord(9,0) = b;
    pointCoord(9,1) = a;
    pointCoord(9,2) = b;

    pointCoord(10,0) = b;
    pointCoord(10,1) = b;
    pointCoord(10,2) = a;

    return pointCoord(i,j);
}

//------------------------------------------------------------------------------
//-------------------------QUADRATURE POINTS - WEIGHTS--------------------------
//------------------------------------------------------------------------------
template<>
double PartQuadrature<2>::WeightList(int i){

pointWeight(	0	) = 	0.00703125000000000 	;
pointWeight(	1	) = 	0.00393559939202584 	;
pointWeight(	2	) = 	0.00393559939202584 	;
pointWeight(	3	) = 	0.00393559939202584 	;
pointWeight(	4	) = 	0.00413731727464081 	;
pointWeight(	5	) = 	0.00413731727464081 	;
pointWeight(	6	) = 	0.00413731727464081 	;
pointWeight(	7	) = 	0.00703125000000000 	;
pointWeight(	8	) = 	0.00393559939202584 	;
pointWeight(	9	) = 	0.00393559939202584 	;
pointWeight(	10	) = 	0.00393559939202584 	;
pointWeight(	11	) = 	0.00413731727464081 	;
pointWeight(	12	) = 	0.00413731727464081 	;
pointWeight(	13	) = 	0.00413731727464081 	;
pointWeight(	14	) = 	0.00703125000000000 	;
pointWeight(	15	) = 	0.00393559939202584 	;
pointWeight(	16	) = 	0.00393559939202584 	;
pointWeight(	17	) = 	0.00393559939202584 	;
pointWeight(	18	) = 	0.00413731727464081 	;
pointWeight(	19	) = 	0.00413731727464081 	;
pointWeight(	20	) = 	0.00413731727464081 	;
pointWeight(	21	) = 	0.00703125000000000 	;
pointWeight(	22	) = 	0.00393559939202584 	;
pointWeight(	23	) = 	0.00393559939202584 	;
pointWeight(	24	) = 	0.00393559939202584 	;
pointWeight(	25	) = 	0.00413731727464081 	;
pointWeight(	26	) = 	0.00413731727464081 	;
pointWeight(	27	) = 	0.00413731727464081 	;
pointWeight(	28	) = 	0.00703125000000000 	;
pointWeight(	29	) = 	0.00393559939202584 	;
pointWeight(	30	) = 	0.00393559939202584 	;
pointWeight(	31	) = 	0.00393559939202584 	;
pointWeight(	32	) = 	0.00413731727464081 	;
pointWeight(	33	) = 	0.00413731727464081 	;
pointWeight(	34	) = 	0.00413731727464081 	;
pointWeight(	35	) = 	0.00703125000000000 	;
pointWeight(	36	) = 	0.00393559939202584 	;
pointWeight(	37	) = 	0.00393559939202584 	;
pointWeight(	38	) = 	0.00393559939202584 	;
pointWeight(	39	) = 	0.00413731727464081 	;
pointWeight(	40	) = 	0.00413731727464081 	;
pointWeight(	41	) = 	0.00413731727464081 	;
pointWeight(	42	) = 	0.00703125000000000 	;
pointWeight(	43	) = 	0.00393559939202584 	;
pointWeight(	44	) = 	0.00393559939202584 	;
pointWeight(	45	) = 	0.00393559939202584 	;
pointWeight(	46	) = 	0.00413731727464081 	;
pointWeight(	47	) = 	0.00413731727464081 	;
pointWeight(	48	) = 	0.00413731727464081 	;
pointWeight(	49	) = 	0.00703125000000000 	;
pointWeight(	50	) = 	0.00393559939202584 	;
pointWeight(	51	) = 	0.00393559939202584 	;
pointWeight(	52	) = 	0.00393559939202584 	;
pointWeight(	53	) = 	0.00413731727464081 	;
pointWeight(	54	) = 	0.00413731727464081 	;
pointWeight(	55	) = 	0.00413731727464081 	;
pointWeight(	56	) = 	0.00703125000000000 	;
pointWeight(	57	) = 	0.00393559939202584 	;
pointWeight(	58	) = 	0.00393559939202584 	;
pointWeight(	59	) = 	0.00393559939202584 	;
pointWeight(	60	) = 	0.00413731727464081 	;
pointWeight(	61	) = 	0.00413731727464081 	;
pointWeight(	62	) = 	0.00413731727464081 	;
pointWeight(	63	) = 	0.00703125000000000 	;
pointWeight(	64	) = 	0.00393559939202584 	;
pointWeight(	65	) = 	0.00393559939202584 	;
pointWeight(	66	) = 	0.00393559939202584 	;
pointWeight(	67	) = 	0.00413731727464081 	;
pointWeight(	68	) = 	0.00413731727464081 	;
pointWeight(	69	) = 	0.00413731727464081 	;
pointWeight(	70	) = 	0.00703125000000000 	;
pointWeight(	71	) = 	0.00393559939202584 	;
pointWeight(	72	) = 	0.00393559939202584 	;
pointWeight(	73	) = 	0.00393559939202584 	;
pointWeight(	74	) = 	0.00413731727464081 	;
pointWeight(	75	) = 	0.00413731727464081 	;
pointWeight(	76	) = 	0.00413731727464081 	;
pointWeight(	77	) = 	0.00703125000000000 	;
pointWeight(	78	) = 	0.00393559939202584 	;
pointWeight(	79	) = 	0.00393559939202584 	;
pointWeight(	80	) = 	0.00393559939202584 	;
pointWeight(	81	) = 	0.00413731727464081 	;
pointWeight(	82	) = 	0.00413731727464081 	;
pointWeight(	83	) = 	0.00413731727464081 	;
pointWeight(	84	) = 	0.00703125000000000 	;
pointWeight(	85	) = 	0.00393559939202584 	;
pointWeight(	86	) = 	0.00393559939202584 	;
pointWeight(	87	) = 	0.00393559939202584 	;
pointWeight(	88	) = 	0.00413731727464081 	;
pointWeight(	89	) = 	0.00413731727464081 	;
pointWeight(	90	) = 	0.00413731727464081 	;
pointWeight(	91	) = 	0.00703125000000000 	;
pointWeight(	92	) = 	0.00393559939202584 	;
pointWeight(	93	) = 	0.00393559939202584 	;
pointWeight(	94	) = 	0.00393559939202584 	;
pointWeight(	95	) = 	0.00413731727464081 	;
pointWeight(	96	) = 	0.00413731727464081 	;
pointWeight(	97	) = 	0.00413731727464081 	;
pointWeight(	98	) = 	0.00703125000000000 	;
pointWeight(	99	) = 	0.00393559939202584 	;
pointWeight(	100	) = 	0.00393559939202584 	;
pointWeight(	101	) = 	0.00393559939202584 	;
pointWeight(	102	) = 	0.00413731727464081 	;
pointWeight(	103	) = 	0.00413731727464081 	;
pointWeight(	104	) = 	0.00413731727464081 	;
pointWeight(	105	) = 	0.00703125000000000 	;
pointWeight(	106	) = 	0.00393559939202584 	;
pointWeight(	107	) = 	0.00393559939202584 	;
pointWeight(	108	) = 	0.00393559939202584 	;
pointWeight(	109	) = 	0.00413731727464081 	;
pointWeight(	110	) = 	0.00413731727464081 	;
pointWeight(	111	) = 	0.00413731727464081 	;

    return pointWeight(i);
};

template<>
double PartQuadrature<3>::WeightList(int i){
    
    pointWeight(0) = -74. / 5625.;
    pointWeight(1) = 343. / 45000.;
    pointWeight(2) = 343. / 45000.;
    pointWeight(3) = 343. / 45000.;
    pointWeight(4) = 343. / 45000.;
    pointWeight(5) = 56. / 2250.;
    pointWeight(6) = 56. / 2250.; 
    pointWeight(7) = 56. / 2250.;
    pointWeight(8) = 56. / 2250.; 
    pointWeight(9) = 56. / 2250.;
    pointWeight(10) = 56. / 2250.; 

    return pointWeight(i);
};

//------------------------------------------------------------------------------
//-----------COMPUTES THE VALUE INTERPOLATED IN THE INTEGRATION POINT-----------
//------------------------------------------------------------------------------
template<>
double PartQuadrature<2>::interpolateQuadraticVariable(NodalValuesQuad nValues,
                                                        int point){
    
    ublas::bounded_vector<double, 2> xsi;

    double int_value = 0.;

    xsi(0) = PointList(point,0);
    xsi(1) = PointList(point,1);
    
    shapeQuad.evaluate(xsi,phi_);
    
    for (int i = 0; i < 6; i++){
        int_value += nValues(i) * phi_(i);
    };

    return int_value;
};

template<>
double PartQuadrature<3>::interpolateQuadraticVariable(NodalValuesQuad nValues,
                                                        int point){
    
    ublas::bounded_vector<double, 3> xsi;

    double int_value = 0.;

    xsi(0) = PointList(point,0);
    xsi(1) = PointList(point,1);
    xsi(2) = PointList(point,2);
    
    shapeQuad.evaluate(xsi,phi_);
    
    for (int i = 0; i < 10; i++){
        int_value += nValues(i) * phi_(i);
    };

    return int_value;
};

//------------------------------------------------------------------------------
//-----------COMPUTES THE VALUE INTERPOLATED IN THE INTEGRATION POINT-----------
//------------------------------------------------------------------------------
template<>
double PartQuadrature<2>::interpolateLinearVariable(NodalValuesLin nValues, 
                                                     int point){
    
    ublas::bounded_vector<double, 2> xsi;

    double int_value = 0.;

    xsi(0) = PointList(point,0);
    xsi(1) = PointList(point,1);
    
    shapeLin.evaluate(xsi,phip_);
    
    for (int i = 0; i < 3; i++){
        int_value += nValues(i) * phip_(i);
    };

    return int_value;
};

#endif
