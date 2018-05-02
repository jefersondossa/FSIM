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
    typedef ublas::bounded_matrix<double, 192,DIM>  PointCoord;

    //Defines the type "PointWeight" which stores the
    //integration points weights
    typedef ublas::bounded_vector<double, 192>      PointWeight;

    //Integration point logical vector
    typedef ublas::bounded_vector<bool, 192>        PointLogical;

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
    pointCoord(	0	,	0	)=	0.0833333333333335	;
pointCoord(	1	,	0	)=	0.0208333333333333	;
pointCoord(	2	,	0	)=	0.0208333333333333	;
pointCoord(	3	,	0	)=	0.2083333333333330	;
pointCoord(	4	,	0	)=	0.1458333333333340	;
pointCoord(	5	,	0	)=	0.1458333333333340	;
pointCoord(	6	,	0	)=	0.0833333333333335	;
pointCoord(	7	,	0	)=	0.0208333333333333	;
pointCoord(	8	,	0	)=	0.0208333333333333	;
pointCoord(	9	,	0	)=	0.0416666666666667	;
pointCoord(	10	,	0	)=	0.1041666666666670	;
pointCoord(	11	,	0	)=	0.1041666666666670	;
pointCoord(	12	,	0	)=	0.3333333333333330	;
pointCoord(	13	,	0	)=	0.2708333333333330	;
pointCoord(	14	,	0	)=	0.2708333333333330	;
pointCoord(	15	,	0	)=	0.4583333333333340	;
pointCoord(	16	,	0	)=	0.3958333333333330	;
pointCoord(	17	,	0	)=	0.3958333333333330	;
pointCoord(	18	,	0	)=	0.3333333333333330	;
pointCoord(	19	,	0	)=	0.2708333333333330	;
pointCoord(	20	,	0	)=	0.2708333333333330	;
pointCoord(	21	,	0	)=	0.2916666666666660	;
pointCoord(	22	,	0	)=	0.3541666666666670	;
pointCoord(	23	,	0	)=	0.3541666666666670	;
pointCoord(	24	,	0	)=	0.0833333333333335	;
pointCoord(	25	,	0	)=	0.0208333333333333	;
pointCoord(	26	,	0	)=	0.0208333333333333	;
pointCoord(	27	,	0	)=	0.2083333333333330	;
pointCoord(	28	,	0	)=	0.1458333333333340	;
pointCoord(	29	,	0	)=	0.1458333333333340	;
pointCoord(	30	,	0	)=	0.0833333333333335	;
pointCoord(	31	,	0	)=	0.0208333333333333	;
pointCoord(	32	,	0	)=	0.0208333333333333	;
pointCoord(	33	,	0	)=	0.0416666666666667	;
pointCoord(	34	,	0	)=	0.1041666666666670	;
pointCoord(	35	,	0	)=	0.1041666666666670	;
pointCoord(	36	,	0	)=	0.1666666666666660	;
pointCoord(	37	,	0	)=	0.2291666666666660	;
pointCoord(	38	,	0	)=	0.2291666666666660	;
pointCoord(	39	,	0	)=	0.0416666666666667	;
pointCoord(	40	,	0	)=	0.1041666666666670	;
pointCoord(	41	,	0	)=	0.1041666666666670	;
pointCoord(	42	,	0	)=	0.1666666666666660	;
pointCoord(	43	,	0	)=	0.2291666666666660	;
pointCoord(	44	,	0	)=	0.2291666666666660	;
pointCoord(	45	,	0	)=	0.2083333333333330	;
pointCoord(	46	,	0	)=	0.1458333333333340	;
pointCoord(	47	,	0	)=	0.1458333333333340	;
pointCoord(	48	,	0	)=	0.5833333333333340	;
pointCoord(	49	,	0	)=	0.5208333333333330	;
pointCoord(	50	,	0	)=	0.5208333333333330	;
pointCoord(	51	,	0	)=	0.7083333333333340	;
pointCoord(	52	,	0	)=	0.6458333333333340	;
pointCoord(	53	,	0	)=	0.6458333333333340	;
pointCoord(	54	,	0	)=	0.5833333333333340	;
pointCoord(	55	,	0	)=	0.5208333333333330	;
pointCoord(	56	,	0	)=	0.5208333333333330	;
pointCoord(	57	,	0	)=	0.5416666666666670	;
pointCoord(	58	,	0	)=	0.6041666666666670	;
pointCoord(	59	,	0	)=	0.6041666666666670	;
pointCoord(	60	,	0	)=	0.8333333333333330	;
pointCoord(	61	,	0	)=	0.7708333333333340	;
pointCoord(	62	,	0	)=	0.7708333333333340	;
pointCoord(	63	,	0	)=	0.9583333333333330	;
pointCoord(	64	,	0	)=	0.8958333333333330	;
pointCoord(	65	,	0	)=	0.8958333333333330	;
pointCoord(	66	,	0	)=	0.8333333333333330	;
pointCoord(	67	,	0	)=	0.7708333333333340	;
pointCoord(	68	,	0	)=	0.7708333333333340	;
pointCoord(	69	,	0	)=	0.7916666666666670	;
pointCoord(	70	,	0	)=	0.8541666666666660	;
pointCoord(	71	,	0	)=	0.8541666666666660	;
pointCoord(	72	,	0	)=	0.5833333333333340	;
pointCoord(	73	,	0	)=	0.5208333333333330	;
pointCoord(	74	,	0	)=	0.5208333333333330	;
pointCoord(	75	,	0	)=	0.7083333333333340	;
pointCoord(	76	,	0	)=	0.6458333333333340	;
pointCoord(	77	,	0	)=	0.6458333333333340	;
pointCoord(	78	,	0	)=	0.5833333333333340	;
pointCoord(	79	,	0	)=	0.5208333333333330	;
pointCoord(	80	,	0	)=	0.5208333333333330	;
pointCoord(	81	,	0	)=	0.5416666666666670	;
pointCoord(	82	,	0	)=	0.6041666666666670	;
pointCoord(	83	,	0	)=	0.6041666666666670	;
pointCoord(	84	,	0	)=	0.6666666666666670	;
pointCoord(	85	,	0	)=	0.7291666666666670	;
pointCoord(	86	,	0	)=	0.7291666666666670	;
pointCoord(	87	,	0	)=	0.5416666666666670	;
pointCoord(	88	,	0	)=	0.6041666666666670	;
pointCoord(	89	,	0	)=	0.6041666666666670	;
pointCoord(	90	,	0	)=	0.6666666666666670	;
pointCoord(	91	,	0	)=	0.7291666666666670	;
pointCoord(	92	,	0	)=	0.7291666666666670	;
pointCoord(	93	,	0	)=	0.7083333333333340	;
pointCoord(	94	,	0	)=	0.6458333333333340	;
pointCoord(	95	,	0	)=	0.6458333333333340	;
pointCoord(	96	,	0	)=	0.4166666666666670	;
pointCoord(	97	,	0	)=	0.4791666666666670	;
pointCoord(	98	,	0	)=	0.4791666666666670	;
pointCoord(	99	,	0	)=	0.2916666666666660	;
pointCoord(	100	,	0	)=	0.3541666666666670	;
pointCoord(	101	,	0	)=	0.3541666666666670	;
pointCoord(	102	,	0	)=	0.4166666666666670	;
pointCoord(	103	,	0	)=	0.4791666666666670	;
pointCoord(	104	,	0	)=	0.4791666666666670	;
pointCoord(	105	,	0	)=	0.4583333333333330	;
pointCoord(	106	,	0	)=	0.3958333333333330	;
pointCoord(	107	,	0	)=	0.3958333333333330	;
pointCoord(	108	,	0	)=	0.1666666666666670	;
pointCoord(	109	,	0	)=	0.2291666666666670	;
pointCoord(	110	,	0	)=	0.2291666666666670	;
pointCoord(	111	,	0	)=	0.0416666666666665	;
pointCoord(	112	,	0	)=	0.1041666666666670	;
pointCoord(	113	,	0	)=	0.1041666666666670	;
pointCoord(	114	,	0	)=	0.1666666666666670	;
pointCoord(	115	,	0	)=	0.2291666666666670	;
pointCoord(	116	,	0	)=	0.2291666666666670	;
pointCoord(	117	,	0	)=	0.2083333333333330	;
pointCoord(	118	,	0	)=	0.1458333333333330	;
pointCoord(	119	,	0	)=	0.1458333333333330	;
pointCoord(	120	,	0	)=	0.4166666666666670	;
pointCoord(	121	,	0	)=	0.4791666666666670	;
pointCoord(	122	,	0	)=	0.4791666666666670	;
pointCoord(	123	,	0	)=	0.2916666666666660	;
pointCoord(	124	,	0	)=	0.3541666666666670	;
pointCoord(	125	,	0	)=	0.3541666666666670	;
pointCoord(	126	,	0	)=	0.4166666666666670	;
pointCoord(	127	,	0	)=	0.4791666666666670	;
pointCoord(	128	,	0	)=	0.4791666666666670	;
pointCoord(	129	,	0	)=	0.4583333333333330	;
pointCoord(	130	,	0	)=	0.3958333333333330	;
pointCoord(	131	,	0	)=	0.3958333333333330	;
pointCoord(	132	,	0	)=	0.3333333333333330	;
pointCoord(	133	,	0	)=	0.2708333333333330	;
pointCoord(	134	,	0	)=	0.2708333333333330	;
pointCoord(	135	,	0	)=	0.4583333333333330	;
pointCoord(	136	,	0	)=	0.3958333333333330	;
pointCoord(	137	,	0	)=	0.3958333333333330	;
pointCoord(	138	,	0	)=	0.3333333333333330	;
pointCoord(	139	,	0	)=	0.2708333333333330	;
pointCoord(	140	,	0	)=	0.2708333333333330	;
pointCoord(	141	,	0	)=	0.2916666666666660	;
pointCoord(	142	,	0	)=	0.3541666666666670	;
pointCoord(	143	,	0	)=	0.3541666666666670	;
pointCoord(	144	,	0	)=	0.0833333333333335	;
pointCoord(	145	,	0	)=	0.0208333333333333	;
pointCoord(	146	,	0	)=	0.0208333333333333	;
pointCoord(	147	,	0	)=	0.2083333333333330	;
pointCoord(	148	,	0	)=	0.1458333333333340	;
pointCoord(	149	,	0	)=	0.1458333333333340	;
pointCoord(	150	,	0	)=	0.0833333333333335	;
pointCoord(	151	,	0	)=	0.0208333333333333	;
pointCoord(	152	,	0	)=	0.0208333333333333	;
pointCoord(	153	,	0	)=	0.0416666666666667	;
pointCoord(	154	,	0	)=	0.1041666666666670	;
pointCoord(	155	,	0	)=	0.1041666666666670	;
pointCoord(	156	,	0	)=	0.3333333333333330	;
pointCoord(	157	,	0	)=	0.2708333333333330	;
pointCoord(	158	,	0	)=	0.2708333333333330	;
pointCoord(	159	,	0	)=	0.4583333333333340	;
pointCoord(	160	,	0	)=	0.3958333333333330	;
pointCoord(	161	,	0	)=	0.3958333333333330	;
pointCoord(	162	,	0	)=	0.3333333333333330	;
pointCoord(	163	,	0	)=	0.2708333333333330	;
pointCoord(	164	,	0	)=	0.2708333333333330	;
pointCoord(	165	,	0	)=	0.2916666666666660	;
pointCoord(	166	,	0	)=	0.3541666666666670	;
pointCoord(	167	,	0	)=	0.3541666666666670	;
pointCoord(	168	,	0	)=	0.0833333333333335	;
pointCoord(	169	,	0	)=	0.0208333333333333	;
pointCoord(	170	,	0	)=	0.0208333333333333	;
pointCoord(	171	,	0	)=	0.2083333333333330	;
pointCoord(	172	,	0	)=	0.1458333333333340	;
pointCoord(	173	,	0	)=	0.1458333333333340	;
pointCoord(	174	,	0	)=	0.0833333333333335	;
pointCoord(	175	,	0	)=	0.0208333333333333	;
pointCoord(	176	,	0	)=	0.0208333333333333	;
pointCoord(	177	,	0	)=	0.0416666666666667	;
pointCoord(	178	,	0	)=	0.1041666666666670	;
pointCoord(	179	,	0	)=	0.1041666666666670	;
pointCoord(	180	,	0	)=	0.1666666666666660	;
pointCoord(	181	,	0	)=	0.2291666666666660	;
pointCoord(	182	,	0	)=	0.2291666666666660	;
pointCoord(	183	,	0	)=	0.0416666666666667	;
pointCoord(	184	,	0	)=	0.1041666666666670	;
pointCoord(	185	,	0	)=	0.1041666666666670	;
pointCoord(	186	,	0	)=	0.1666666666666660	;
pointCoord(	187	,	0	)=	0.2291666666666660	;
pointCoord(	188	,	0	)=	0.2291666666666660	;
pointCoord(	189	,	0	)=	0.2083333333333330	;
pointCoord(	190	,	0	)=	0.1458333333333340	;
pointCoord(	191	,	0	)=	0.1458333333333340	;
pointCoord(	0	,	1	)=	0.0208333333333333	;
pointCoord(	1	,	1	)=	0.0833333333333335	;
pointCoord(	2	,	1	)=	0.0208333333333333	;
pointCoord(	3	,	1	)=	0.0208333333333333	;
pointCoord(	4	,	1	)=	0.0833333333333335	;
pointCoord(	5	,	1	)=	0.0208333333333333	;
pointCoord(	6	,	1	)=	0.1458333333333340	;
pointCoord(	7	,	1	)=	0.2083333333333330	;
pointCoord(	8	,	1	)=	0.1458333333333340	;
pointCoord(	9	,	1	)=	0.1041666666666670	;
pointCoord(	10	,	1	)=	0.0416666666666667	;
pointCoord(	11	,	1	)=	0.1041666666666670	;
pointCoord(	12	,	1	)=	0.0208333333333333	;
pointCoord(	13	,	1	)=	0.0833333333333335	;
pointCoord(	14	,	1	)=	0.0208333333333333	;
pointCoord(	15	,	1	)=	0.0208333333333333	;
pointCoord(	16	,	1	)=	0.0833333333333335	;
pointCoord(	17	,	1	)=	0.0208333333333333	;
pointCoord(	18	,	1	)=	0.1458333333333340	;
pointCoord(	19	,	1	)=	0.2083333333333330	;
pointCoord(	20	,	1	)=	0.1458333333333340	;
pointCoord(	21	,	1	)=	0.1041666666666670	;
pointCoord(	22	,	1	)=	0.0416666666666667	;
pointCoord(	23	,	1	)=	0.1041666666666670	;
pointCoord(	24	,	1	)=	0.2708333333333330	;
pointCoord(	25	,	1	)=	0.3333333333333330	;
pointCoord(	26	,	1	)=	0.2708333333333330	;
pointCoord(	27	,	1	)=	0.2708333333333330	;
pointCoord(	28	,	1	)=	0.3333333333333330	;
pointCoord(	29	,	1	)=	0.2708333333333330	;
pointCoord(	30	,	1	)=	0.3958333333333330	;
pointCoord(	31	,	1	)=	0.4583333333333340	;
pointCoord(	32	,	1	)=	0.3958333333333330	;
pointCoord(	33	,	1	)=	0.3541666666666670	;
pointCoord(	34	,	1	)=	0.2916666666666660	;
pointCoord(	35	,	1	)=	0.3541666666666670	;
pointCoord(	36	,	1	)=	0.2291666666666660	;
pointCoord(	37	,	1	)=	0.1666666666666660	;
pointCoord(	38	,	1	)=	0.2291666666666660	;
pointCoord(	39	,	1	)=	0.2291666666666660	;
pointCoord(	40	,	1	)=	0.1666666666666660	;
pointCoord(	41	,	1	)=	0.2291666666666660	;
pointCoord(	42	,	1	)=	0.1041666666666670	;
pointCoord(	43	,	1	)=	0.0416666666666667	;
pointCoord(	44	,	1	)=	0.1041666666666670	;
pointCoord(	45	,	1	)=	0.1458333333333340	;
pointCoord(	46	,	1	)=	0.2083333333333330	;
pointCoord(	47	,	1	)=	0.1458333333333340	;
pointCoord(	48	,	1	)=	0.0208333333333333	;
pointCoord(	49	,	1	)=	0.0833333333333335	;
pointCoord(	50	,	1	)=	0.0208333333333333	;
pointCoord(	51	,	1	)=	0.0208333333333333	;
pointCoord(	52	,	1	)=	0.0833333333333335	;
pointCoord(	53	,	1	)=	0.0208333333333333	;
pointCoord(	54	,	1	)=	0.1458333333333340	;
pointCoord(	55	,	1	)=	0.2083333333333330	;
pointCoord(	56	,	1	)=	0.1458333333333340	;
pointCoord(	57	,	1	)=	0.1041666666666670	;
pointCoord(	58	,	1	)=	0.0416666666666667	;
pointCoord(	59	,	1	)=	0.1041666666666670	;
pointCoord(	60	,	1	)=	0.0208333333333333	;
pointCoord(	61	,	1	)=	0.0833333333333335	;
pointCoord(	62	,	1	)=	0.0208333333333333	;
pointCoord(	63	,	1	)=	0.0208333333333333	;
pointCoord(	64	,	1	)=	0.0833333333333335	;
pointCoord(	65	,	1	)=	0.0208333333333333	;
pointCoord(	66	,	1	)=	0.1458333333333340	;
pointCoord(	67	,	1	)=	0.2083333333333330	;
pointCoord(	68	,	1	)=	0.1458333333333340	;
pointCoord(	69	,	1	)=	0.1041666666666670	;
pointCoord(	70	,	1	)=	0.0416666666666667	;
pointCoord(	71	,	1	)=	0.1041666666666670	;
pointCoord(	72	,	1	)=	0.2708333333333330	;
pointCoord(	73	,	1	)=	0.3333333333333330	;
pointCoord(	74	,	1	)=	0.2708333333333330	;
pointCoord(	75	,	1	)=	0.2708333333333330	;
pointCoord(	76	,	1	)=	0.3333333333333330	;
pointCoord(	77	,	1	)=	0.2708333333333330	;
pointCoord(	78	,	1	)=	0.3958333333333330	;
pointCoord(	79	,	1	)=	0.4583333333333340	;
pointCoord(	80	,	1	)=	0.3958333333333330	;
pointCoord(	81	,	1	)=	0.3541666666666670	;
pointCoord(	82	,	1	)=	0.2916666666666660	;
pointCoord(	83	,	1	)=	0.3541666666666670	;
pointCoord(	84	,	1	)=	0.2291666666666660	;
pointCoord(	85	,	1	)=	0.1666666666666660	;
pointCoord(	86	,	1	)=	0.2291666666666660	;
pointCoord(	87	,	1	)=	0.2291666666666660	;
pointCoord(	88	,	1	)=	0.1666666666666660	;
pointCoord(	89	,	1	)=	0.2291666666666660	;
pointCoord(	90	,	1	)=	0.1041666666666670	;
pointCoord(	91	,	1	)=	0.0416666666666667	;
pointCoord(	92	,	1	)=	0.1041666666666670	;
pointCoord(	93	,	1	)=	0.1458333333333340	;
pointCoord(	94	,	1	)=	0.2083333333333330	;
pointCoord(	95	,	1	)=	0.1458333333333340	;
pointCoord(	96	,	1	)=	0.4791666666666670	;
pointCoord(	97	,	1	)=	0.4166666666666670	;
pointCoord(	98	,	1	)=	0.4791666666666670	;
pointCoord(	99	,	1	)=	0.4791666666666670	;
pointCoord(	100	,	1	)=	0.4166666666666670	;
pointCoord(	101	,	1	)=	0.4791666666666670	;
pointCoord(	102	,	1	)=	0.3541666666666670	;
pointCoord(	103	,	1	)=	0.2916666666666660	;
pointCoord(	104	,	1	)=	0.3541666666666670	;
pointCoord(	105	,	1	)=	0.3958333333333330	;
pointCoord(	106	,	1	)=	0.4583333333333330	;
pointCoord(	107	,	1	)=	0.3958333333333330	;
pointCoord(	108	,	1	)=	0.4791666666666670	;
pointCoord(	109	,	1	)=	0.4166666666666670	;
pointCoord(	110	,	1	)=	0.4791666666666670	;
pointCoord(	111	,	1	)=	0.4791666666666670	;
pointCoord(	112	,	1	)=	0.4166666666666670	;
pointCoord(	113	,	1	)=	0.4791666666666670	;
pointCoord(	114	,	1	)=	0.3541666666666670	;
pointCoord(	115	,	1	)=	0.2916666666666660	;
pointCoord(	116	,	1	)=	0.3541666666666670	;
pointCoord(	117	,	1	)=	0.3958333333333330	;
pointCoord(	118	,	1	)=	0.4583333333333330	;
pointCoord(	119	,	1	)=	0.3958333333333330	;
pointCoord(	120	,	1	)=	0.2291666666666670	;
pointCoord(	121	,	1	)=	0.1666666666666670	;
pointCoord(	122	,	1	)=	0.2291666666666670	;
pointCoord(	123	,	1	)=	0.2291666666666670	;
pointCoord(	124	,	1	)=	0.1666666666666670	;
pointCoord(	125	,	1	)=	0.2291666666666670	;
pointCoord(	126	,	1	)=	0.1041666666666670	;
pointCoord(	127	,	1	)=	0.0416666666666665	;
pointCoord(	128	,	1	)=	0.1041666666666670	;
pointCoord(	129	,	1	)=	0.1458333333333330	;
pointCoord(	130	,	1	)=	0.2083333333333330	;
pointCoord(	131	,	1	)=	0.1458333333333330	;
pointCoord(	132	,	1	)=	0.2708333333333330	;
pointCoord(	133	,	1	)=	0.3333333333333330	;
pointCoord(	134	,	1	)=	0.2708333333333330	;
pointCoord(	135	,	1	)=	0.2708333333333330	;
pointCoord(	136	,	1	)=	0.3333333333333330	;
pointCoord(	137	,	1	)=	0.2708333333333330	;
pointCoord(	138	,	1	)=	0.3958333333333330	;
pointCoord(	139	,	1	)=	0.4583333333333330	;
pointCoord(	140	,	1	)=	0.3958333333333330	;
pointCoord(	141	,	1	)=	0.3541666666666670	;
pointCoord(	142	,	1	)=	0.2916666666666660	;
pointCoord(	143	,	1	)=	0.3541666666666670	;
pointCoord(	144	,	1	)=	0.5208333333333330	;
pointCoord(	145	,	1	)=	0.5833333333333340	;
pointCoord(	146	,	1	)=	0.5208333333333330	;
pointCoord(	147	,	1	)=	0.5208333333333330	;
pointCoord(	148	,	1	)=	0.5833333333333340	;
pointCoord(	149	,	1	)=	0.5208333333333330	;
pointCoord(	150	,	1	)=	0.6458333333333340	;
pointCoord(	151	,	1	)=	0.7083333333333340	;
pointCoord(	152	,	1	)=	0.6458333333333340	;
pointCoord(	153	,	1	)=	0.6041666666666670	;
pointCoord(	154	,	1	)=	0.5416666666666670	;
pointCoord(	155	,	1	)=	0.6041666666666670	;
pointCoord(	156	,	1	)=	0.5208333333333330	;
pointCoord(	157	,	1	)=	0.5833333333333340	;
pointCoord(	158	,	1	)=	0.5208333333333330	;
pointCoord(	159	,	1	)=	0.5208333333333330	;
pointCoord(	160	,	1	)=	0.5833333333333340	;
pointCoord(	161	,	1	)=	0.5208333333333330	;
pointCoord(	162	,	1	)=	0.6458333333333340	;
pointCoord(	163	,	1	)=	0.7083333333333340	;
pointCoord(	164	,	1	)=	0.6458333333333340	;
pointCoord(	165	,	1	)=	0.6041666666666670	;
pointCoord(	166	,	1	)=	0.5416666666666670	;
pointCoord(	167	,	1	)=	0.6041666666666670	;
pointCoord(	168	,	1	)=	0.7708333333333340	;
pointCoord(	169	,	1	)=	0.8333333333333330	;
pointCoord(	170	,	1	)=	0.7708333333333340	;
pointCoord(	171	,	1	)=	0.7708333333333340	;
pointCoord(	172	,	1	)=	0.8333333333333330	;
pointCoord(	173	,	1	)=	0.7708333333333340	;
pointCoord(	174	,	1	)=	0.8958333333333330	;
pointCoord(	175	,	1	)=	0.9583333333333330	;
pointCoord(	176	,	1	)=	0.8958333333333330	;
pointCoord(	177	,	1	)=	0.8541666666666660	;
pointCoord(	178	,	1	)=	0.7916666666666670	;
pointCoord(	179	,	1	)=	0.8541666666666660	;
pointCoord(	180	,	1	)=	0.7291666666666670	;
pointCoord(	181	,	1	)=	0.6666666666666670	;
pointCoord(	182	,	1	)=	0.7291666666666670	;
pointCoord(	183	,	1	)=	0.7291666666666670	;
pointCoord(	184	,	1	)=	0.6666666666666670	;
pointCoord(	185	,	1	)=	0.7291666666666670	;
pointCoord(	186	,	1	)=	0.6041666666666670	;
pointCoord(	187	,	1	)=	0.5416666666666670	;
pointCoord(	188	,	1	)=	0.6041666666666670	;
pointCoord(	189	,	1	)=	0.6458333333333340	;
pointCoord(	190	,	1	)=	0.7083333333333340	;
pointCoord(	191	,	1	)=	0.6458333333333340	;
 
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
    pointWeight(	0	)=	1./384.	;
pointWeight(	1	)=	1./384.	;
pointWeight(	2	)=	1./384.	;
pointWeight(	3	)=	1./384.	;
pointWeight(	4	)=	1./384.	;
pointWeight(	5	)=	1./384.	;
pointWeight(	6	)=	1./384.	;
pointWeight(	7	)=	1./384.	;
pointWeight(	8	)=	1./384.	;
pointWeight(	9	)=	1./384.	;
pointWeight(	10	)=	1./384.	;
pointWeight(	11	)=	1./384.	;
pointWeight(	12	)=	1./384.	;
pointWeight(	13	)=	1./384.	;
pointWeight(	14	)=	1./384.	;
pointWeight(	15	)=	1./384.	;
pointWeight(	16	)=	1./384.	;
pointWeight(	17	)=	1./384.	;
pointWeight(	18	)=	1./384.	;
pointWeight(	19	)=	1./384.	;
pointWeight(	20	)=	1./384.	;
pointWeight(	21	)=	1./384.	;
pointWeight(	22	)=	1./384.	;
pointWeight(	23	)=	1./384.	;
pointWeight(	24	)=	1./384.	;
pointWeight(	25	)=	1./384.	;
pointWeight(	26	)=	1./384.	;
pointWeight(	27	)=	1./384.	;
pointWeight(	28	)=	1./384.	;
pointWeight(	29	)=	1./384.	;
pointWeight(	30	)=	1./384.	;
pointWeight(	31	)=	1./384.	;
pointWeight(	32	)=	1./384.	;
pointWeight(	33	)=	1./384.	;
pointWeight(	34	)=	1./384.	;
pointWeight(	35	)=	1./384.	;
pointWeight(	36	)=	1./384.	;
pointWeight(	37	)=	1./384.	;
pointWeight(	38	)=	1./384.	;
pointWeight(	39	)=	1./384.	;
pointWeight(	40	)=	1./384.	;
pointWeight(	41	)=	1./384.	;
pointWeight(	42	)=	1./384.	;
pointWeight(	43	)=	1./384.	;
pointWeight(	44	)=	1./384.	;
pointWeight(	45	)=	1./384.	;
pointWeight(	46	)=	1./384.	;
pointWeight(	47	)=	1./384.	;
pointWeight(	48	)=	1./384.	;
pointWeight(	49	)=	1./384.	;
pointWeight(	50	)=	1./384.	;
pointWeight(	51	)=	1./384.	;
pointWeight(	52	)=	1./384.	;
pointWeight(	53	)=	1./384.	;
pointWeight(	54	)=	1./384.	;
pointWeight(	55	)=	1./384.	;
pointWeight(	56	)=	1./384.	;
pointWeight(	57	)=	1./384.	;
pointWeight(	58	)=	1./384.	;
pointWeight(	59	)=	1./384.	;
pointWeight(	60	)=	1./384.	;
pointWeight(	61	)=	1./384.	;
pointWeight(	62	)=	1./384.	;
pointWeight(	63	)=	1./384.	;
pointWeight(	64	)=	1./384.	;
pointWeight(	65	)=	1./384.	;
pointWeight(	66	)=	1./384.	;
pointWeight(	67	)=	1./384.	;
pointWeight(	68	)=	1./384.	;
pointWeight(	69	)=	1./384.	;
pointWeight(	70	)=	1./384.	;
pointWeight(	71	)=	1./384.	;
pointWeight(	72	)=	1./384.	;
pointWeight(	73	)=	1./384.	;
pointWeight(	74	)=	1./384.	;
pointWeight(	75	)=	1./384.	;
pointWeight(	76	)=	1./384.	;
pointWeight(	77	)=	1./384.	;
pointWeight(	78	)=	1./384.	;
pointWeight(	79	)=	1./384.	;
pointWeight(	80	)=	1./384.	;
pointWeight(	81	)=	1./384.	;
pointWeight(	82	)=	1./384.	;
pointWeight(	83	)=	1./384.	;
pointWeight(	84	)=	1./384.	;
pointWeight(	85	)=	1./384.	;
pointWeight(	86	)=	1./384.	;
pointWeight(	87	)=	1./384.	;
pointWeight(	88	)=	1./384.	;
pointWeight(	89	)=	1./384.	;
pointWeight(	90	)=	1./384.	;
pointWeight(	91	)=	1./384.	;
pointWeight(	92	)=	1./384.	;
pointWeight(	93	)=	1./384.	;
pointWeight(	94	)=	1./384.	;
pointWeight(	95	)=	1./384.	;
pointWeight(	96	)=	1./384.	;
pointWeight(	97	)=	1./384.	;
pointWeight(	98	)=	1./384.	;
pointWeight(	99	)=	1./384.	;
pointWeight(	100	)=	1./384.	;
pointWeight(	101	)=	1./384.	;
pointWeight(	102	)=	1./384.	;
pointWeight(	103	)=	1./384.	;
pointWeight(	104	)=	1./384.	;
pointWeight(	105	)=	1./384.	;
pointWeight(	106	)=	1./384.	;
pointWeight(	107	)=	1./384.	;
pointWeight(	108	)=	1./384.	;
pointWeight(	109	)=	1./384.	;
pointWeight(	110	)=	1./384.	;
pointWeight(	111	)=	1./384.	;
pointWeight(	112	)=	1./384.	;
pointWeight(	113	)=	1./384.	;
pointWeight(	114	)=	1./384.	;
pointWeight(	115	)=	1./384.	;
pointWeight(	116	)=	1./384.	;
pointWeight(	117	)=	1./384.	;
pointWeight(	118	)=	1./384.	;
pointWeight(	119	)=	1./384.	;
pointWeight(	120	)=	1./384.	;
pointWeight(	121	)=	1./384.	;
pointWeight(	122	)=	1./384.	;
pointWeight(	123	)=	1./384.	;
pointWeight(	124	)=	1./384.	;
pointWeight(	125	)=	1./384.	;
pointWeight(	126	)=	1./384.	;
pointWeight(	127	)=	1./384.	;
pointWeight(	128	)=	1./384.	;
pointWeight(	129	)=	1./384.	;
pointWeight(	130	)=	1./384.	;
pointWeight(	131	)=	1./384.	;
pointWeight(	132	)=	1./384.	;
pointWeight(	133	)=	1./384.	;
pointWeight(	134	)=	1./384.	;
pointWeight(	135	)=	1./384.	;
pointWeight(	136	)=	1./384.	;
pointWeight(	137	)=	1./384.	;
pointWeight(	138	)=	1./384.	;
pointWeight(	139	)=	1./384.	;
pointWeight(	140	)=	1./384.	;
pointWeight(	141	)=	1./384.	;
pointWeight(	142	)=	1./384.	;
pointWeight(	143	)=	1./384.	;
pointWeight(	144	)=	1./384.	;
pointWeight(	145	)=	1./384.	;
pointWeight(	146	)=	1./384.	;
pointWeight(	147	)=	1./384.	;
pointWeight(	148	)=	1./384.	;
pointWeight(	149	)=	1./384.	;
pointWeight(	150	)=	1./384.	;
pointWeight(	151	)=	1./384.	;
pointWeight(	152	)=	1./384.	;
pointWeight(	153	)=	1./384.	;
pointWeight(	154	)=	1./384.	;
pointWeight(	155	)=	1./384.	;
pointWeight(	156	)=	1./384.	;
pointWeight(	157	)=	1./384.	;
pointWeight(	158	)=	1./384.	;
pointWeight(	159	)=	1./384.	;
pointWeight(	160	)=	1./384.	;
pointWeight(	161	)=	1./384.	;
pointWeight(	162	)=	1./384.	;
pointWeight(	163	)=	1./384.	;
pointWeight(	164	)=	1./384.	;
pointWeight(	165	)=	1./384.	;
pointWeight(	166	)=	1./384.	;
pointWeight(	167	)=	1./384.	;
pointWeight(	168	)=	1./384.	;
pointWeight(	169	)=	1./384.	;
pointWeight(	170	)=	1./384.	;
pointWeight(	171	)=	1./384.	;
pointWeight(	172	)=	1./384.	;
pointWeight(	173	)=	1./384.	;
pointWeight(	174	)=	1./384.	;
pointWeight(	175	)=	1./384.	;
pointWeight(	176	)=	1./384.	;
pointWeight(	177	)=	1./384.	;
pointWeight(	178	)=	1./384.	;
pointWeight(	179	)=	1./384.	;
pointWeight(	180	)=	1./384.	;
pointWeight(	181	)=	1./384.	;
pointWeight(	182	)=	1./384.	;
pointWeight(	183	)=	1./384.	;
pointWeight(	184	)=	1./384.	;
pointWeight(	185	)=	1./384.	;
pointWeight(	186	)=	1./384.	;
pointWeight(	187	)=	1./384.	;
pointWeight(	188	)=	1./384.	;
pointWeight(	189	)=	1./384.	;
pointWeight(	190	)=	1./384.	;
pointWeight(	191	)=	1./384.	;

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
