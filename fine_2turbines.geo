p0 = newp; Point(p0) = {0.000000, 0.000000, 0.000000, 0.100000};
//
p1 = newp; Point(p1) = {0.000000, 1.850000, 0.000000, 0.100000};
//
p2 = newp; Point(p2) = {0.000000, -1.850000, 0.000000, 0.100000};
//
p3 = newp; Point(p3) = {0.000000, -3.850000, 0.000000, 0.100000};
//
p4 = newp; Point(p4) = {0.000000, 0.150000, 0.000000, 0.100000};
//
p5 = newp; Point(p5) = {0.000000, -0.150000, 0.000000, 0.100000};
//
p6 = newp; Point(p6) = {0.000000, 3.850000, 0.000000, 0.100000};
//
l0 = newl; Circle(l0) = {p4, p2, p3}; Physical Line('l0') = {l0};
//
l1 = newl; Circle(l1) = {p5, p1, p6}; Physical Line('l1') = {l1};
//
l2 = newl; Circle(l2) = {p4, p2, p3}; Physical Line('l2') = {l2};
//
l3 = newl; Circle(l3) = {p5, p1, p6}; Physical Line('l3') = {l3};
//
l4 = newl; Circle(l4) = {p4, p0, p5}; Physical Line('l4') = {l4};
//
l5 = newl; Circle(l5) = {p5, p0, p4}; Physical Line('l5') = {l5};
//
p7 = newp; Point(p7) = {-0.955065, 5.416443, 0.000000, 0.100000};
//
p8 = newp; Point(p8) = {-3.939231, -0.694593, 0.000000, 0.100000};
//
p9 = newp; Point(p9) = {0.955065, -5.416443, 0.000000, 0.100000};
//
p10 = newp; Point(p10) = {3.939231, 0.694593, 0.000000, 0.100000};
//
l6 = newl; Line(l6) = {p6, p7}; Physical Line('l6') = {l6};
//
l7 = newl; Line(l7) = {p3, p9}; Physical Line('l7') = {l7};
//
p11 = newp; Point(p11) = {-1.128713, 6.401250, 0.000000, 0.100000};
//
p12 = newp; Point(p12) = {-4.924039, -0.868241, 0.000000, 0.100000};
//
p13 = newp; Point(p13) = {1.128713, -6.401250, 0.000000, 0.100000};
//
p14 = newp; Point(p14) = {4.924039, 0.868241, 0.000000, 0.100000};
//
l8 = newl; Ellipse(l8) = {p7, p0, p7, p8}; Physical Line('l8') = {l8};
//
l9 = newl; Ellipse(l9) = {p8, p0, p8, p9}; Physical Line('l9') = {l9};
//
l10 = newl; Ellipse(l10) = {p9, p0, p9, p10}; Physical Line('l10') = {l10};
//
l11 = newl; Ellipse(l11) = {p10, p0, p10, p7}; Physical Line('l11') = {l11};
//
l12 = newl; Ellipse(l12) = {p11, p0, p11, p12}; Physical Line('l12') = {l12};
//
l13 = newl; Ellipse(l13) = {p12, p0, p12, p13}; Physical Line('l13') = {l13};
//
l14 = newl; Ellipse(l14) = {p13, p0, p13, p14}; Physical Line('l14') = {l14};
//
l15 = newl; Ellipse(l15) = {p14, p0, p14, p11}; Physical Line('l15') = {l15};
//
Transfinite Line {l8, l9, l10, l11, l12, l13, l14, l15} = 36 Using Progression 1;
//
ll0 = newll; Line Loop(ll0) = {l8, l9, -l7, -l0, -l5, l1, l6};
//
ll1 = newll; Line Loop(ll1) = {-l8, -l9, -l10, -l11, l12, l13, l14, l15};
//
ll2 = newll; Line Loop(ll2) = {l10, l11, -l6, -l3, -l4, l2, l7};
//
s0 = news; Plane Surface(s0) = {ll0}; Physical Surface('s0') = {s0};
//
s1 = news; Plane Surface(s1) = {ll1}; Physical Surface('s1') = {s1};
//
s2 = news; Plane Surface(s2) = {ll2}; Physical Surface('s2') = {s2};
//
p15 = newp; Point(p15) = {14.000000, 0.000000, 0.000000, 0.100000};
//
p16 = newp; Point(p16) = {14.000000, 1.850000, 0.000000, 0.100000};
//
p17 = newp; Point(p17) = {14.000000, -1.850000, 0.000000, 0.100000};
//
p18 = newp; Point(p18) = {14.000000, -3.850000, 0.000000, 0.100000};
//
p19 = newp; Point(p19) = {14.000000, 0.150000, 0.000000, 0.100000};
//
p20 = newp; Point(p20) = {14.000000, -0.150000, 0.000000, 0.100000};
//
p21 = newp; Point(p21) = {14.000000, 3.850000, 0.000000, 0.100000};
//
l16 = newl; Circle(l16) = {p19, p17, p18}; Physical Line('l16') = {l16};
//
l17 = newl; Circle(l17) = {p20, p16, p21}; Physical Line('l17') = {l17};
//
l18 = newl; Circle(l18) = {p19, p17, p18}; Physical Line('l18') = {l18};
//
l19 = newl; Circle(l19) = {p20, p16, p21}; Physical Line('l19') = {l19};
//
l20 = newl; Circle(l20) = {p19, p15, p20}; Physical Line('l20') = {l20};
//
l21 = newl; Circle(l21) = {p20, p15, p19}; Physical Line('l21') = {l21};
//
p22 = newp; Point(p22) = {13.044935, 5.416443, 0.000000, 0.100000};
//
p23 = newp; Point(p23) = {10.060769, -0.694593, 0.000000, 0.100000};
//
p24 = newp; Point(p24) = {14.955065, -5.416443, 0.000000, 0.100000};
//
p25 = newp; Point(p25) = {17.939231, 0.694593, 0.000000, 0.100000};
//
l22 = newl; Line(l22) = {p21, p22}; Physical Line('l22') = {l22};
//
l23 = newl; Line(l23) = {p18, p24}; Physical Line('l23') = {l23};
//
p26 = newp; Point(p26) = {12.871287, 6.401250, 0.000000, 0.100000};
//
p27 = newp; Point(p27) = {9.075961, -0.868241, 0.000000, 0.100000};
//
p28 = newp; Point(p28) = {15.128713, -6.401250, 0.000000, 0.100000};
//
p29 = newp; Point(p29) = {18.924039, 0.868241, 0.000000, 0.100000};
//
l24 = newl; Ellipse(l24) = {p22, p15, p22, p23}; Physical Line('l24') = {l24};
//
l25 = newl; Ellipse(l25) = {p23, p15, p23, p24}; Physical Line('l25') = {l25};
//
l26 = newl; Ellipse(l26) = {p24, p15, p24, p25}; Physical Line('l26') = {l26};
//
l27 = newl; Ellipse(l27) = {p25, p15, p25, p22}; Physical Line('l27') = {l27};
//
l28 = newl; Ellipse(l28) = {p26, p15, p26, p27}; Physical Line('l28') = {l28};
//
l29 = newl; Ellipse(l29) = {p27, p15, p27, p28}; Physical Line('l29') = {l29};
//
l30 = newl; Ellipse(l30) = {p28, p15, p28, p29}; Physical Line('l30') = {l30};
//
l31 = newl; Ellipse(l31) = {p29, p15, p29, p26}; Physical Line('l31') = {l31};
//
Transfinite Line {l24, l25, l26, l27, l28, l29, l30, l31} = 36 Using Progression 1;
//
ll3 = newll; Line Loop(ll3) = {l24, l25, -l23, -l16, -l21, l17, l22};
//
ll4 = newll; Line Loop(ll4) = {-l24, -l25, -l26, -l27, l28, l29, l30, l31};
//
ll5 = newll; Line Loop(ll5) = {l26, l27, -l22, -l19, -l20, l18, l23};
//
s3 = news; Plane Surface(s3) = {ll3}; Physical Surface('s3') = {s3};
//
s4 = news; Plane Surface(s4) = {ll4}; Physical Surface('s4') = {s4};
//
s5 = news; Plane Surface(s5) = {ll5}; Physical Surface('s5') = {s5};
//
