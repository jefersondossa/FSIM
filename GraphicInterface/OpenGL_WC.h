#ifndef OPENGL_WC_H
#define OPENGL_WC_H

#include <FL/Fl.H>
#include <FL/Fl_Gl_Window.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Group.H>
#include <FL/Fl_Menu_Bar.H>
#include <FL/gl.h>
#include <FL/glu.h>
#include <FL/glut.H>
#include <FL/fl_message.H>
#include <iostream>
#include <fstream>
#include <vector>
#include <tuple>
#include <regex>

using namespace std;


void IdleCallback(void* pData);
class Playback : public Fl_Gl_Window {

private:
    int last_x, last_y;
    bool rotating;
    float rotationX, rotationY;
    float zoom;

    //vector<tuple<int, float, float, float>> points;
    //vector<tuple<int, int>> lines;

    vector <float> points;
    vector <float> lines;
    vector <float> circles;
    string geopath;
    string mshpath;
    vector<int> matids;


    double sum_x;
    double sum_y;
    double sum_z;

    double average_x;
    double average_y;
    double average_z;

public:

    void draw() override;
    void InitializeGL();
    //void ReadGeometry(const string filename);
    void DrawGeometry();
    void DrawAxes();
    virtual int handle(int event) override;

    float rotation;
    float rotationIncrement;
    bool  animating;

public:
    Playback(int X, int Y, int W, int H, const char *L);

    void SetGeoPath(string geop){
        geopath = geop;
    }
    void SetMshPath(string mshp){
        mshpath = mshp;
    }

    void Setmatid(vector<int> matID){

        matids = matID;
    }
};
#endif
