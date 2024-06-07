/* This cpp file describes the following Weak Form Menu Inline and Static functions:

2. From fApply_WF_cb:
        fInline_Apply_WF;
        fStatic_Apply_WF;

    From fElasticity2D_cb:
        fInline_Elasticity2D;
        fStatic_Elasticity2D;

    From fElasticityPositional2D_cb:
        fInline_ElasticityPositional2D;
        fStatic_ElasticityPositional2D;

    From fElasticity3D_cb:
        fInline_Elasticity3D;
        fStatic_Elasticity3D;

    From fElasticTruss_cb:
        fInline_ElasticTruss;
        fStatic_ElasticTruss;

    From fPositionalTruss_cb:
        fInline_PositionalTruss;
        fStatic_PositionalTruss;

    From fStokes_cb:
        fInline_Stokes;
        fStatic_Stokes;

    From fNavierStokes_cb:
        fInline_NavierStokes;
        fStatic_NavierStokes;

    From fPoisson_cb:
        fInline_Poisson;
        fStatic_Poisson;

    From fL2Projection_cb:
        fInline_L2Projection;
        fStatic_L2Projection;

*/

#include "WindowConstructor.h"

/* ========================= From fApply_WF_cb ========================= */

void WindowConstructor::fInline_Apply_WF(Fl_Button*, void*){
    this->fApply_WF_cb();
}

void WindowConstructor::fStatic_Apply_WF(Fl_Button* o, void* v){
    ((WindowConstructor*)(o->parent()->parent()->user_data()))->fInline_Apply_WF(o,v);
}

/* ========================= From fElasticity2D_cb ========================= */

void WindowConstructor::fInline_Elasticity2D(Fl_Menu_*, void*){
    this->fElasticity2D_cb();
}

void WindowConstructor::fStatic_Elasticity2D(Fl_Menu_* o, void* v){
    ((WindowConstructor*)(o->parent()->parent()->user_data()))->fInline_Elasticity2D(o,v);
}

/* ========================= From fElasticityPositional2D_cb ========================= */

void WindowConstructor::fInline_ElasticityPositional2D(Fl_Menu_*, void*){
    this->fElasticityPositional2D_cb();
}

void WindowConstructor::fStatic_ElasticityPositional2D(Fl_Menu_* o, void* v){
    ((WindowConstructor*)(o->parent()->parent()->user_data()))->fInline_ElasticityPositional2D(o,v);
}

/* ========================= From fElasticity3D_cb ========================= */

void WindowConstructor::fInline_Elasticity3D(Fl_Menu_*, void*){
    this->fElasticity3D_cb();
}

void WindowConstructor::fStatic_Elasticity3D(Fl_Menu_* o, void* v){
    ((WindowConstructor*)(o->parent()->parent()->user_data()))->fInline_Elasticity3D(o,v);
}

/* ========================= From fElasticTruss_cb ========================= */

void WindowConstructor::fInline_ElasticTruss(Fl_Menu_*, void*){
    this->fElasticTruss_cb();
}

void WindowConstructor::fStatic_ElasticTruss(Fl_Menu_* o, void* v){
    ((WindowConstructor*)(o->parent()->parent()->user_data()))->fInline_ElasticTruss(o,v);
}

/* ========================= From fPositionalTruss_cb ========================= */

void WindowConstructor::fInline_PositionalTruss(Fl_Menu_*, void*){
    this->fPositionalTruss_cb();
}

void WindowConstructor::fStatic_PositionalTruss(Fl_Menu_* o, void* v){
    ((WindowConstructor*)(o->parent()->parent()->user_data()))->fInline_PositionalTruss(o,v);
}

/* ========================= From fStokes_cb ========================= */

void WindowConstructor::fInline_Stokes(Fl_Menu_*, void*){
    this->fStokes_cb();
}

void WindowConstructor::fStatic_Stokes(Fl_Menu_* o, void* v){
    ((WindowConstructor*)(o->parent()->parent()->user_data()))->fInline_Stokes(o,v);
}


/* ========================= From fNavierStokes_cb ========================= */
void WindowConstructor::fInline_NavierStokes(Fl_Menu_*, void*){
    this->fNavierStokes_cb();
}

void WindowConstructor::fStatic_NavierStokes(Fl_Menu_* o, void* v){
    ((WindowConstructor*)(o->parent()->parent()->user_data()))->fInline_NavierStokes(o,v);
}

/* ========================= From fPoisson_cb ========================= */

void WindowConstructor::fInline_Poisson(Fl_Menu_*, void*){
    this->fPoisson_cb();
}

void WindowConstructor::fStatic_Poisson(Fl_Menu_* o, void* v){
    ((WindowConstructor*)(o->parent()->parent()->user_data()))->fInline_Poisson(o,v);
}

/* ========================= From fL2Projection_cb ========================= */

void WindowConstructor::fInline_L2Projection(Fl_Menu_*, void*){
    this->fL2Projection_cb();
}

void WindowConstructor::fStatic_L2Projection(Fl_Menu_* o, void* v){
    ((WindowConstructor*)(o->parent()->parent()->user_data()))->fInline_L2Projection(o,v);
}
