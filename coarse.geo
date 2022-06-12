p0 = newp; Point(p0) = {-50.000000, -60.000000, 0.000000, 2.000000};
//
p1 = newp; Point(p1) = {-9.000000, -60.000000, 0.000000, 2.000000};
//
p2 = newp; Point(p2) = {9.000000, -60.000000, 0.000000, 2.000000};
//
p3 = newp; Point(p3) = {250.000000, -60.000000, 0.000000, 2.000000};
//
p4 = newp; Point(p4) = {-50.000000, -9.000000, 0.000000, 2.000000};
//
p5 = newp; Point(p5) = {-9.000000, -9.000000, 0.000000, 0.133333};
//
p6 = newp; Point(p6) = {9.000000, -9.000000, 0.000000, 0.133333};
//
p7 = newp; Point(p7) = {250.000000, -9.000000, 0.000000, 2.000000};
//
p8 = newp; Point(p8) = {-50.000000, 9.000000, 0.000000, 2.000000};
//
p9 = newp; Point(p9) = {-9.000000, 9.000000, 0.000000, 0.133333};
//
p10 = newp; Point(p10) = {9.000000, 9.000000, 0.000000, 0.133333};
//
p11 = newp; Point(p11) = {250.000000, 9.000000, 0.000000, 2.000000};
//
p12 = newp; Point(p12) = {-50.000000, 60.000000, 0.000000, 2.000000};
//
p13 = newp; Point(p13) = {-9.000000, 60.000000, 0.000000, 2.000000};
//
p14 = newp; Point(p14) = {9.000000, 60.000000, 0.000000, 2.000000};
//
p15 = newp; Point(p15) = {250.000000, 60.000000, 0.000000, 2.000000};
//
l0 = newl; Line(l0) = {p0, p1}; Physical Line('l0') = {l0};
//
l1 = newl; Line(l1) = {p1, p2}; Physical Line('l1') = {l1};
//
l2 = newl; Line(l2) = {p2, p3}; Physical Line('l2') = {l2};
//
l3 = newl; Line(l3) = {p0, p4}; Physical Line('l3') = {l3};
//
l4 = newl; Line(l4) = {p1, p5}; Physical Line('l4') = {l4};
//
l5 = newl; Line(l5) = {p2, p6}; Physical Line('l5') = {l5};
//
l6 = newl; Line(l6) = {p3, p7}; Physical Line('l6') = {l6};
//
l7 = newl; Line(l7) = {p4, p5}; Physical Line('l7') = {l7};
//
l8 = newl; Line(l8) = {p5, p6}; Physical Line('l8') = {l8};
//
l9 = newl; Line(l9) = {p6, p7}; Physical Line('l9') = {l9};
//
l10 = newl; Line(l10) = {p4, p8}; Physical Line('l10') = {l10};
//
l11 = newl; Line(l11) = {p5, p9}; Physical Line('l11') = {l11};
//
l12 = newl; Line(l12) = {p6, p10}; Physical Line('l12') = {l12};
//
l13 = newl; Line(l13) = {p7, p11}; Physical Line('l13') = {l13};
//
l14 = newl; Line(l14) = {p8, p9}; Physical Line('l14') = {l14};
//
l15 = newl; Line(l15) = {p9, p10}; Physical Line('l15') = {l15};
//
l16 = newl; Line(l16) = {p10, p11}; Physical Line('l16') = {l16};
//
l17 = newl; Line(l17) = {p8, p12}; Physical Line('l17') = {l17};
//
l18 = newl; Line(l18) = {p9, p13}; Physical Line('l18') = {l18};
//
l19 = newl; Line(l19) = {p10, p14}; Physical Line('l19') = {l19};
//
l20 = newl; Line(l20) = {p11, p15}; Physical Line('l20') = {l20};
//
l21 = newl; Line(l21) = {p12, p13}; Physical Line('l21') = {l21};
//
l22 = newl; Line(l22) = {p13, p14}; Physical Line('l22') = {l22};
//
l23 = newl; Line(l23) = {p14, p15}; Physical Line('l23') = {l23};
//
ll0 = newll; Line Loop(ll0) = {l0, l4, -l7, -l3};
//
ll1 = newll; Line Loop(ll1) = {l1, l5, -l8, -l4};
//
ll2 = newll; Line Loop(ll2) = {l2, l6, -l9, -l5};
//
ll3 = newll; Line Loop(ll3) = {l7, l11, -l14, -l10};
//
ll4 = newll; Line Loop(ll4) = {l8, l12, -l15, -l11};
//
ll5 = newll; Line Loop(ll5) = {l9, l13, -l16, -l12};
//
ll6 = newll; Line Loop(ll6) = {l14, l18, -l21, -l17};
//
ll7 = newll; Line Loop(ll7) = {l15, l19, -l22, -l18};
//
ll8 = newll; Line Loop(ll8) = {l16, l20, -l23, -l19};
//
s0 = news; Plane Surface(s0) = {ll0}; Physical Surface('s0') = {s0};
//
s1 = news; Plane Surface(s1) = {ll1}; Physical Surface('s1') = {s1};
//
s2 = news; Plane Surface(s2) = {ll2}; Physical Surface('s2') = {s2};
//
s3 = news; Plane Surface(s3) = {ll3}; Physical Surface('s3') = {s3};
//
s4 = news; Plane Surface(s4) = {ll4}; Physical Surface('s4') = {s4};
//
s5 = news; Plane Surface(s5) = {ll5}; Physical Surface('s5') = {s5};
//
s6 = news; Plane Surface(s6) = {ll6}; Physical Surface('s6') = {s6};
//
s7 = news; Plane Surface(s7) = {ll7}; Physical Surface('s7') = {s7};
//
s8 = news; Plane Surface(s8) = {ll8}; Physical Surface('s8') = {s8};
//
