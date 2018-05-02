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

    pointCoord(	0	,	0	) =	0.16666666666666700	;
    pointCoord(	1	,	0	) =	0.47401085907171200	;
    pointCoord(	2	,	0	) =	0.01299457046414420	;
    pointCoord(	3	,	0	) =	0.01299457046414420	;
    pointCoord(	4	,	0	) =	0.40571249735207700	;
    pointCoord(	5	,	0	) =	0.04714375132396140	;
    pointCoord(	6	,	0	) =	0.04714375132396140	;
    pointCoord(	7	,	0	) =	0.00536322498278530	;
    pointCoord(	8	,	0	) =	0.24731838750860700	;
    pointCoord(	9	,	0	) =	0.24731838750860700	;
    pointCoord(	10	,	0	) =	0.29265661738548600	;
    pointCoord(	11	,	0	) =	0.10367169130725700	;
    pointCoord(	12	,	0	) =	0.10367169130725700	;
    pointCoord(	13	,	0	) =	0.06109219429950940	;
    pointCoord(	14	,	0	) =	0.21945390285024500	;
    pointCoord(	15	,	0	) =	0.21945390285024500	;
    pointCoord(	16	,	0	) =	0.33896882744129500	;
    pointCoord(	17	,	0	) =	0.33896882744129500	;
    pointCoord(	18	,	0	) =	0.02242083879456530	;
    pointCoord(	19	,	0	) =	0.02242083879456530	;
    pointCoord(	20	,	0	) =	0.13861033376414000	;
    pointCoord(	21	,	0	) =	0.13861033376414000	;
    pointCoord(	22	,	0	) =	0.42943514064131800	;
    pointCoord(	23	,	0	) =	0.42943514064131800	;
    pointCoord(	24	,	0	) =	0.00000000000000000	;
    pointCoord(	25	,	0	) =	0.00000000000000000	;
    pointCoord(	26	,	0	) =	0.07056485935868180	;
    pointCoord(	27	,	0	) =	0.07056485935868180	;
    pointCoord(	28	,	0	) =	0.66666666666666700	;
    pointCoord(	29	,	0	) =	0.97401085907171200	;
    pointCoord(	30	,	0	) =	0.51299457046414400	;
    pointCoord(	31	,	0	) =	0.51299457046414400	;
    pointCoord(	32	,	0	) =	0.90571249735207700	;
    pointCoord(	33	,	0	) =	0.54714375132396100	;
    pointCoord(	34	,	0	) =	0.54714375132396100	;
    pointCoord(	35	,	0	) =	0.50536322498278500	;
    pointCoord(	36	,	0	) =	0.74731838750860700	;
    pointCoord(	37	,	0	) =	0.74731838750860700	;
    pointCoord(	38	,	0	) =	0.79265661738548600	;
    pointCoord(	39	,	0	) =	0.60367169130725700	;
    pointCoord(	40	,	0	) =	0.60367169130725700	;
    pointCoord(	41	,	0	) =	0.56109219429950900	;
    pointCoord(	42	,	0	) =	0.71945390285024500	;
    pointCoord(	43	,	0	) =	0.71945390285024500	;
    pointCoord(	44	,	0	) =	0.83896882744129500	;
    pointCoord(	45	,	0	) =	0.83896882744129500	;
    pointCoord(	46	,	0	) =	0.52242083879456500	;
    pointCoord(	47	,	0	) =	0.52242083879456500	;
    pointCoord(	48	,	0	) =	0.63861033376414000	;
    pointCoord(	49	,	0	) =	0.63861033376414000	;
    pointCoord(	50	,	0	) =	0.92943514064131800	;
    pointCoord(	51	,	0	) =	0.92943514064131800	;
    pointCoord(	52	,	0	) =	0.50000000000000000	;
    pointCoord(	53	,	0	) =	0.50000000000000000	;
    pointCoord(	54	,	0	) =	0.57056485935868200	;
    pointCoord(	55	,	0	) =	0.57056485935868200	;
    pointCoord(	56	,	0	) =	0.16666666666666700	;
    pointCoord(	57	,	0	) =	0.47401085907171200	;
    pointCoord(	58	,	0	) =	0.01299457046414420	;
    pointCoord(	59	,	0	) =	0.01299457046414420	;
    pointCoord(	60	,	0	) =	0.40571249735207700	;
    pointCoord(	61	,	0	) =	0.04714375132396140	;
    pointCoord(	62	,	0	) =	0.04714375132396140	;
    pointCoord(	63	,	0	) =	0.00536322498278530	;
    pointCoord(	64	,	0	) =	0.24731838750860700	;
    pointCoord(	65	,	0	) =	0.24731838750860700	;
pointCoord(	66	,	0	) =	0.29265661738548600	;
pointCoord(	67	,	0	) =	0.10367169130725700	;
pointCoord(	68	,	0	) =	0.10367169130725700	;
pointCoord(	69	,	0	) =	0.06109219429950940	;
pointCoord(	70	,	0	) =	0.21945390285024500	;
pointCoord(	71	,	0	) =	0.21945390285024500	;
pointCoord(	72	,	0	) =	0.33896882744129500	;
pointCoord(	73	,	0	) =	0.33896882744129500	;
pointCoord(	74	,	0	) =	0.02242083879456530	;
pointCoord(	75	,	0	) =	0.02242083879456530	;
pointCoord(	76	,	0	) =	0.13861033376414000	;
pointCoord(	77	,	0	) =	0.13861033376414000	;
pointCoord(	78	,	0	) =	0.42943514064131800	;
pointCoord(	79	,	0	) =	0.42943514064131800	;
pointCoord(	80	,	0	) =	0.00000000000000000	;
pointCoord(	81	,	0	) =	0.00000000000000000	;
pointCoord(	82	,	0	) =	0.07056485935868180	;
pointCoord(	83	,	0	) =	0.07056485935868180	;
pointCoord(	84	,	0	) =	0.33333333333333300	;
pointCoord(	85	,	0	) =	0.02598914092828840	;
pointCoord(	86	,	0	) =	0.48700542953585600	;
pointCoord(	87	,	0	) =	0.48700542953585600	;
pointCoord(	88	,	0	) =	0.09428750264792270	;
pointCoord(	89	,	0	) =	0.45285624867603900	;
pointCoord(	90	,	0	) =	0.45285624867603900	;
pointCoord(	91	,	0	) =	0.49463677501721500	;
pointCoord(	92	,	0	) =	0.25268161249139300	;
pointCoord(	93	,	0	) =	0.25268161249139300	;
pointCoord(	94	,	0	) =	0.20734338261451400	;
pointCoord(	95	,	0	) =	0.39632830869274300	;
pointCoord(	96	,	0	) =	0.39632830869274300	;
pointCoord(	97	,	0	) =	0.43890780570049100	;
pointCoord(	98	,	0	) =	0.28054609714975500	;
pointCoord(	99	,	0	) =	0.28054609714975500	;
pointCoord(	100	,	0	) =	0.16103117255870500	;
pointCoord(	101	,	0	) =	0.16103117255870500	;
pointCoord(	102	,	0	) =	0.47757916120543500	;
pointCoord(	103	,	0	) =	0.47757916120543500	;
pointCoord(	104	,	0	) =	0.36138966623586000	;
pointCoord(	105	,	0	) =	0.36138966623586000	;
pointCoord(	106	,	0	) =	0.07056485935868180	;
pointCoord(	107	,	0	) =	0.07056485935868180	;
pointCoord(	108	,	0	) =	0.50000000000000000	;
pointCoord(	109	,	0	) =	0.50000000000000000	;
pointCoord(	110	,	0	) =	0.42943514064131800	;
pointCoord(	111	,	0	) =	0.42943514064131800	;
pointCoord(	0	,	1	) =	0.16666666666666700	;
pointCoord(	1	,	1	) =	0.01299457046414420	;
pointCoord(	2	,	1	) =	0.47401085907171200	;
pointCoord(	3	,	1	) =	0.01299457046414420	;
pointCoord(	4	,	1	) =	0.04714375132396140	;
pointCoord(	5	,	1	) =	0.40571249735207700	;
pointCoord(	6	,	1	) =	0.04714375132396140	;
pointCoord(	7	,	1	) =	0.24731838750860700	;
pointCoord(	8	,	1	) =	0.00536322498278530	;
pointCoord(	9	,	1	) =	0.24731838750860700	;
pointCoord(	10	,	1	) =	0.10367169130725700	;
pointCoord(	11	,	1	) =	0.29265661738548600	;
pointCoord(	12	,	1	) =	0.10367169130725700	;
pointCoord(	13	,	1	) =	0.21945390285024500	;
pointCoord(	14	,	1	) =	0.06109219429950940	;
pointCoord(	15	,	1	) =	0.21945390285024500	;
pointCoord(	16	,	1	) =	0.02242083879456530	;
pointCoord(	17	,	1	) =	0.13861033376414000	;
pointCoord(	18	,	1	) =	0.33896882744129500	;
pointCoord(	19	,	1	) =	0.13861033376414000	;
pointCoord(	20	,	1	) =	0.33896882744129500	;
pointCoord(	21	,	1	) =	0.02242083879456530	;
pointCoord(	22	,	1	) =	0.00000000000000000	;
pointCoord(	23	,	1	) =	0.07056485935868180	;
pointCoord(	24	,	1	) =	0.42943514064131800	;
pointCoord(	25	,	1	) =	0.07056485935868180	;
pointCoord(	26	,	1	) =	0.42943514064131800	;
pointCoord(	27	,	1	) =	0.00000000000000000	;
pointCoord(	28	,	1	) =	0.16666666666666700	;
pointCoord(	29	,	1	) =	0.01299457046414420	;
pointCoord(	30	,	1	) =	0.47401085907171200	;
pointCoord(	31	,	1	) =	0.01299457046414420	;
pointCoord(	32	,	1	) =	0.04714375132396140	;
pointCoord(	33	,	1	) =	0.40571249735207700	;
pointCoord(	34	,	1	) =	0.04714375132396140	;
pointCoord(	35	,	1	) =	0.24731838750860700	;
pointCoord(	36	,	1	) =	0.00536322498278530	;
pointCoord(	37	,	1	) =	0.24731838750860700	;
pointCoord(	38	,	1	) =	0.10367169130725700	;
pointCoord(	39	,	1	) =	0.29265661738548600	;
pointCoord(	40	,	1	) =	0.10367169130725700	;
pointCoord(	41	,	1	) =	0.21945390285024500	;
pointCoord(	42	,	1	) =	0.06109219429950940	;
pointCoord(	43	,	1	) =	0.21945390285024500	;
pointCoord(	44	,	1	) =	0.02242083879456530	;
pointCoord(	45	,	1	) =	0.13861033376414000	;
pointCoord(	46	,	1	) =	0.33896882744129500	;
pointCoord(	47	,	1	) =	0.13861033376414000	;
pointCoord(	48	,	1	) =	0.33896882744129500	;
pointCoord(	49	,	1	) =	0.02242083879456530	;
pointCoord(	50	,	1	) =	0.00000000000000000	;
pointCoord(	51	,	1	) =	0.07056485935868180	;
pointCoord(	52	,	1	) =	0.42943514064131800	;
pointCoord(	53	,	1	) =	0.07056485935868180	;
pointCoord(	54	,	1	) =	0.42943514064131800	;
pointCoord(	55	,	1	) =	0.00000000000000000	;
pointCoord(	56	,	1	) =	0.66666666666666700	;
pointCoord(	57	,	1	) =	0.51299457046414400	;
pointCoord(	58	,	1	) =	0.97401085907171200	;
pointCoord(	59	,	1	) =	0.51299457046414400	;
pointCoord(	60	,	1	) =	0.54714375132396100	;
pointCoord(	61	,	1	) =	0.90571249735207700	;
pointCoord(	62	,	1	) =	0.54714375132396100	;
pointCoord(	63	,	1	) =	0.74731838750860700	;
pointCoord(	64	,	1	) =	0.50536322498278500	;
pointCoord(	65	,	1	) =	0.74731838750860700	;
pointCoord(	66	,	1	) =	0.60367169130725700	;
pointCoord(	67	,	1	) =	0.79265661738548600	;
pointCoord(	68	,	1	) =	0.60367169130725700	;
pointCoord(	69	,	1	) =	0.71945390285024500	;
pointCoord(	70	,	1	) =	0.56109219429950900	;
pointCoord(	71	,	1	) =	0.71945390285024500	;
pointCoord(	72	,	1	) =	0.52242083879456500	;
pointCoord(	73	,	1	) =	0.63861033376414000	;
pointCoord(	74	,	1	) =	0.83896882744129500	;
pointCoord(	75	,	1	) =	0.63861033376414000	;
pointCoord(	76	,	1	) =	0.83896882744129500	;
pointCoord(	77	,	1	) =	0.52242083879456500	;
pointCoord(	78	,	1	) =	0.50000000000000000	;
pointCoord(	79	,	1	) =	0.57056485935868200	;
pointCoord(	80	,	1	) =	0.92943514064131800	;
pointCoord(	81	,	1	) =	0.57056485935868200	;
pointCoord(	82	,	1	) =	0.92943514064131800	;
pointCoord(	83	,	1	) =	0.50000000000000000	;
pointCoord(	84	,	1	) =	0.33333333333333300	;
pointCoord(	85	,	1	) =	0.48700542953585600	;
pointCoord(	86	,	1	) =	0.02598914092828840	;
pointCoord(	87	,	1	) =	0.48700542953585600	;
pointCoord(	88	,	1	) =	0.45285624867603900	;
pointCoord(	89	,	1	) =	0.09428750264792260	;
pointCoord(	90	,	1	) =	0.45285624867603900	;
pointCoord(	91	,	1	) =	0.25268161249139300	;
pointCoord(	92	,	1	) =	0.49463677501721500	;
pointCoord(	93	,	1	) =	0.25268161249139300	;
pointCoord(	94	,	1	) =	0.39632830869274300	;
pointCoord(	95	,	1	) =	0.20734338261451400	;
pointCoord(	96	,	1	) =	0.39632830869274300	;
pointCoord(	97	,	1	) =	0.28054609714975500	;
pointCoord(	98	,	1	) =	0.43890780570049100	;
pointCoord(	99	,	1	) =	0.28054609714975500	;
pointCoord(	100	,	1	) =	0.47757916120543500	;
pointCoord(	101	,	1	) =	0.36138966623586000	;
pointCoord(	102	,	1	) =	0.16103117255870500	;
pointCoord(	103	,	1	) =	0.36138966623586000	;
pointCoord(	104	,	1	) =	0.16103117255870500	;
pointCoord(	105	,	1	) =	0.47757916120543500	;
pointCoord(	106	,	1	) =	0.50000000000000000	;
pointCoord(	107	,	1	) =	0.42943514064131800	;
pointCoord(	108	,	1	) =	0.07056485935868180	;
pointCoord(	109	,	1	) =	0.42943514064131800	;
pointCoord(	110	,	1	) =	0.07056485935868180	;
pointCoord(	111	,	1	) =	0.50000000000000000	;


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
    
    // pointWeight(0) = 0.11250;
    // pointWeight(1) = 0.125939180544827 / 2.;
    // pointWeight(2) = 0.125939180544827 / 2.;
    // pointWeight(3) = 0.125939180544827 / 2.;
    // pointWeight(4) = 0.132394152788506 / 2.;
    // pointWeight(5) = 0.132394152788506 / 2.;
    // pointWeight(6) = 0.132394152788506 / 2.;
    pointWeight(	0	) =	0.010997162645277700	;
pointWeight(	1	) =	0.001093038944217020	;
pointWeight(	2	) =	0.001093038944217020	;
pointWeight(	3	) =	0.001093038944217020	;
pointWeight(	4	) =	0.004760196499241920	;
pointWeight(	5	) =	0.004760196499241920	;
pointWeight(	6	) =	0.004760196499241920	;
pointWeight(	7	) =	0.002356931007016410	;
pointWeight(	8	) =	0.002356931007016410	;
pointWeight(	9	) =	0.002356931007016410	;
pointWeight(	10	) =	0.009019962193092630	;
pointWeight(	11	) =	0.009019962193092630	;
pointWeight(	12	) =	0.009019962193092630	;
pointWeight(	13	) =	0.008666142338192150	;
pointWeight(	14	) =	0.008666142338192150	;
pointWeight(	15	) =	0.008666142338192150	;
pointWeight(	16	) =	0.005132039428661080	;
pointWeight(	17	) =	0.005132039428661080	;
pointWeight(	18	) =	0.005132039428661080	;
pointWeight(	19	) =	0.005132039428661080	;
pointWeight(	20	) =	0.005132039428661080	;
pointWeight(	21	) =	0.005132039428661080	;
pointWeight(	22	) =	0.000920297972912572	;
pointWeight(	23	) =	0.000920297972912572	;
pointWeight(	24	) =	0.000920297972912572	;
pointWeight(	25	) =	0.000920297972912572	;
pointWeight(	26	) =	0.000920297972912572	;
pointWeight(	27	) =	0.000920297972912572	;
pointWeight(	28	) =	0.010997162645277700	;
pointWeight(	29	) =	0.001093038944217020	;
pointWeight(	30	) =	0.001093038944217020	;
pointWeight(	31	) =	0.001093038944217020	;
pointWeight(	32	) =	0.004760196499241920	;
pointWeight(	33	) =	0.004760196499241920	;
pointWeight(	34	) =	0.004760196499241920	;
pointWeight(	35	) =	0.002356931007016410	;
pointWeight(	36	) =	0.002356931007016410	;
pointWeight(	37	) =	0.002356931007016410	;
pointWeight(	38	) =	0.009019962193092630	;
pointWeight(	39	) =	0.009019962193092630	;
pointWeight(	40	) =	0.009019962193092630	;
pointWeight(	41	) =	0.008666142338192150	;
pointWeight(	42	) =	0.008666142338192150	;
pointWeight(	43	) =	0.008666142338192150	;
pointWeight(	44	) =	0.005132039428661080	;
pointWeight(	45	) =	0.005132039428661080	;
pointWeight(	46	) =	0.005132039428661080	;
pointWeight(	47	) =	0.005132039428661080	;
pointWeight(	48	) =	0.005132039428661080	;
pointWeight(	49	) =	0.005132039428661080	;
pointWeight(	50	) =	0.000920297972912572	;
pointWeight(	51	) =	0.000920297972912572	;
pointWeight(	52	) =	0.000920297972912572	;
pointWeight(	53	) =	0.000920297972912572	;
pointWeight(	54	) =	0.000920297972912572	;
pointWeight(	55	) =	0.000920297972912572	;
pointWeight(	56	) =	0.010997162645277700	;
pointWeight(	57	) =	0.001093038944217020	;
pointWeight(	58	) =	0.001093038944217020	;
pointWeight(	59	) =	0.001093038944217020	;
pointWeight(	60	) =	0.004760196499241920	;
pointWeight(	61	) =	0.004760196499241920	;
pointWeight(	62	) =	0.004760196499241920	;
pointWeight(	63	) =	0.002356931007016410	;
pointWeight(	64	) =	0.002356931007016410	;
pointWeight(	65	) =	0.002356931007016410	;
pointWeight(	66	) =	0.009019962193092630	;
pointWeight(	67	) =	0.009019962193092630	;
pointWeight(	68	) =	0.009019962193092630	;
pointWeight(	69	) =	0.008666142338192150	;
pointWeight(	70	) =	0.008666142338192150	;
pointWeight(	71	) =	0.008666142338192150	;
pointWeight(	72	) =	0.005132039428661080	;
pointWeight(	73	) =	0.005132039428661080	;
pointWeight(	74	) =	0.005132039428661080	;
pointWeight(	75	) =	0.005132039428661080	;
pointWeight(	76	) =	0.005132039428661080	;
pointWeight(	77	) =	0.005132039428661080	;
pointWeight(	78	) =	0.000920297972912572	;
pointWeight(	79	) =	0.000920297972912572	;
pointWeight(	80	) =	0.000920297972912572	;
pointWeight(	81	) =	0.000920297972912572	;
pointWeight(	82	) =	0.000920297972912572	;
pointWeight(	83	) =	0.000920297972912572	;
pointWeight(	84	) =	0.010997162645277700	;
pointWeight(	85	) =	0.001093038944217020	;
pointWeight(	86	) =	0.001093038944217020	;
pointWeight(	87	) =	0.001093038944217020	;
pointWeight(	88	) =	0.004760196499241920	;
pointWeight(	89	) =	0.004760196499241920	;
pointWeight(	90	) =	0.004760196499241920	;
pointWeight(	91	) =	0.002356931007016410	;
pointWeight(	92	) =	0.002356931007016410	;
pointWeight(	93	) =	0.002356931007016410	;
pointWeight(	94	) =	0.009019962193092630	;
pointWeight(	95	) =	0.009019962193092630	;
pointWeight(	96	) =	0.009019962193092630	;
pointWeight(	97	) =	0.008666142338192150	;
pointWeight(	98	) =	0.008666142338192150	;
pointWeight(	99	) =	0.008666142338192150	;
pointWeight(	100	) =	0.005132039428661080	;
pointWeight(	101	) =	0.005132039428661080	;
pointWeight(	102	) =	0.005132039428661080	;
pointWeight(	103	) =	0.005132039428661080	;
pointWeight(	104	) =	0.005132039428661080	;
pointWeight(	105	) =	0.005132039428661080	;
pointWeight(	106	) =	0.000920297972912572	;
pointWeight(	107	) =	0.000920297972912572	;
pointWeight(	108	) =	0.000920297972912572	;
pointWeight(	109	) =	0.000920297972912572	;
pointWeight(	110	) =	0.000920297972912572	;
pointWeight(	111	) =	0.000920297972912572	;

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
double PartQuadrature<2>::interpolateQuadraticVariable(
                                                      NodalValuesQuad nValues,
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
double PartQuadrature<3>::interpolateQuadraticVariable(
                                                       NodalValuesQuad nValues,
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
double PartQuadrature<2>::interpolateLinearVariable(
                                                     NodalValuesLin nValues, 
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
