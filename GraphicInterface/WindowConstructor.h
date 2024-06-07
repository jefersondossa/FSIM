/*  This WindowConstructor.h header file contains all objects used to create the FSArl GUI. To achieve this, it uses the FLTK and OpenGL libraries (C++ GUI tools).
    
    This file is subdivided in:
    1. Library includes;
    2. WindowConstructor class:
      2.1. GLOBAL VARIABLES -> A few variables that are used in general;
      2.2. SCRIPT DISPLAY COMPONENTS -> Allows the user to see and verify his choices along the problem creation;
      2.3. GEOMETRY MENU COMPONENTS -> Allows the user to inseart the .geo and .msh files. They are essencial to load the geometry and mesh;
      2.4. WEAK FORM MENU COMPONENTS -> Allows the user to set the material variables and the external actions;
      2.5. ANALYSIS MENU COMPONENTS -> Allows the user to set the Analysis type, whether instantaneous or in steps;
      2.6. POS-PROCESSING MENU COMPONENTS -> Allows the user to set and visualise properly the vector fields and deformations;
      2.7. TOP HEADER MENU COMPONENTS -> The FSArl window header;
      2.8. TOOL BAR MENU COMPONENTS -> Just buttons to access the 2.3 to 2.6 menus; 
    
    Each WindowConstructor class part is also subvided in 3 other parts:
    a. WIDGETS -> The visual components from FSArls window;
    b. CREATION FUNCTIONS -> Defines the widgets configurations;
    c. CALLBACK FUNCTIONS -> Defines the callback from each widget;
*/

#ifndef WindowConstructor_h
#define WindowConstructor_h

//**************************************
//FLTK library and OpenGL window Header.
//**************************************
#include <FL/Fl.H>
#include <FL/Fl_Double_Window.H>
#include <FL/Fl_Group.H>
#include <FL/Fl_Output.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Value_Output.H>
#include <FL/Fl_Choice.H>
#include <FL/Fl_Value_Input.H>
#include <FL/Fl_Check_Button.H>
#include <FL/Fl_Box.H>
#include <FL/Fl_Text_Display.H>
#include <FL/Fl_Input.H>
#include <FL/Fl_Menu_Bar.H>
#include <FL/Fl_File_Chooser.H>
#include <FL/Fl_Image.H>
#include <Fl/Fl_PNG_Image.H>

#include <OpenGLWindow/OpenGL_WC.h>

//*****************
//Default libraries
//*****************
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <regex>
#include <map>


using namespace std;

class WindowConstructor {

//******************************************
//GLOBAL VARIABLES - FSArl GUI
//******************************************

public:

  Fl_Text_Buffer Buffer; //GUI -> Defines the Buffer that stores the Script text;
  fstream ScriptMemory_txt; //GUI -> Defines a txt editor file that haves Script Memory information;
  string Str;
  ostringstream Oss;
  /* Since the variables are contained in objects and, since Fl_Text_Buffer function only accepts "text", it's necessary to save the entire
  text ("text" << objects) inside the Oss object. */

  string Line;
  int Dimension; //GUI -> Saves the problem dimension.
  smatch Match; //GUI -> Declare the "match" name for each pattern found in the desired files:
  string gp; //GUI -> Saves the choosen .geo file path;
  string mp; //GUI -> Saves the choosen .msh file path;
  vector<int> matids; //GUI -> Saves the Phisical Groups id in to a vector;
  vector<string> names; //GUI -> Saves the Physical Groups names in to a vector;
 
  WindowConstructor();  //GUI -> WindowCOnstructor constructor;
  void f_Open(); //GUI -> Defines the functions that opens the Window;
  Fl_Double_Window *Window; //GUI -> Defines the FLTK window;
  Playback *playback; //GUI -> Defines the OpenGL window.


//******************************************
//SCRIPT DISPLAY COMPONENTS - FSArl GUI
//******************************************

  /* ========================= WIDGETS ========================= */

public:

  Fl_Group *Script_Menu; //GUI -> Defines the window that contains all Script Display components;
  Fl_Text_Display *ScriptDisplay; //GUI -> Defines the ScriptDisplay window;
  Fl_Button *Button_SaveScript; //GUI -> Defines the button that opens the Script_Popup_Menu;

  Fl_Group *Script_Popup_Menu; //GUI -> Defines the SaveScript window;
  Fl_Box *Script_Popup_Text;
  Fl_Input *Script_Text_Name; //GUI -> Defines the Input for the txt file that stores the ScriptDisplay text for latter modifications;
  Fl_Button *Button_TextCancel; //GUI -> Defines the button that cancels the save Script Menu action;
  Fl_Button *Button_TextSave; //GUI ->  Defines the button that allows the user to save the Script Display in a txt file.

  /* ========================= CREATION FUNCTIONS ========================= */

public:

  void fScript_Menu(Fl_Double_Window *window); //GUI -> Defines the Script Display window function;
  void fScript_Popup_Menu(Fl_Double_Window *window); //GUI -> Defines the Save Script Display window function.

  /* ========================= CALLBACK FUNCTIONS ========================= */

public:

  void fButton_SaveScript_cb(); //GUI -> Defines the Save Script Button callback;
  void fButton_TextSave_cb(); //GUI -> Defines the Save Script Text button callback;
  void fButton_TextCancel_cb(); //GUI -> Defines the Cancel Script Text button callback.

private:

  inline void fInline_Button_SaveScript(Fl_Button*, void*);
  static void fStatic_Button_SaveScript(Fl_Button*, void*);

  inline void fInline_Button_TextSave(Fl_Button*, void*);
  static void fStatic_Button_TextSave(Fl_Button*, void*);

  inline void fInline_Button_TextCancel(Fl_Button*, void*);
  static void fStatic_Button_TextCancel(Fl_Button*, void*);


//************************************
//GEOMETRY MENU COMPONENTS - FSArl GUI
//************************************

  /* ========================= WIDGETS ========================= */

public:

  Fl_Group *Geometry_Menu; //GUI -> Defines the Geometry Menu;

  Fl_Output *geo; //GUI -> Defines the  geo_path output from Geometry Menu;
  Fl_Output *msh; //GUI -> Defines the  msh_path output (Geometry Menu widget);
  Fl_Value_Output *dimension; //GUI -> Defines the  problem dimension output (Geometry Menu widget);
  Fl_Button *Apply_geometry; //GUI -> Defines the  Save button for the Geometry Menu widget;

  map<string,int> PhysGroup; //GUI -> Allows to list the Physical Groups Nmaes with the Material ID.

  /* ========================= CREATION FUNCTIONS ========================= */

public:

  void fGeometry(Fl_Double_Window *Window); //GUI -> Defines the Geometry Menu function;
  void fOpen_Main_Window(); //GUI -> Defines the function that opens the FSArl window.

  /* ========================= CALLBACK FUNCTIONS ========================= */

public:

  void fWF_free_cb(); //GUI -> Defines the function that frees the Weak Form Submenus utilization. This reduces the chances of the user making mistakes;
  void fInput_geo_cb(); //GUI -> Defines the functionality from geo Fl_Output;
  void fInput_msh_cb(); //GUI -> Defines the functionality from msh Fl_Output;
  void fApply_geometry_cb(); //GUI -> Defines the function that saves the geo and msh paths. It also reads and writes the paths in the script.

private:

  inline void fInline_WF_free(Fl_Output*, void*);
  static void fStatic_WF_free(Fl_Output*, void*);

  inline void fInline_Input_geo(Fl_Output*, void*);
  static void fStatic_Input_geo(Fl_Output*, void*);

  inline void fInline_Input_msh(Fl_Output*, void*);
  static void fStatic_Input_msh(Fl_Output*, void*);

  inline void fInline_Apply_geometry(Fl_Button*, void*);
  static void fStatic_Apply_geometry(Fl_Button*, void*);


//*************************************
//WEAK FORM MENU COMPONENTS - FSArl GUI
//*************************************

  /* ========================= WIDGETS ========================= */

public:

  Fl_Group *WF_Menu; //GUI -> Defines the Weak Form Menu;
  Fl_Group *Elasticity2D_Menu; //GUI -> Defines the Elasticity 2D Submenu from Weak Form Menu;
  Fl_Group *ElasticityPositional2D_Menu; //GUI -> Defines the Elasticity Positional 2D Submenu from Weak Form Menu;
  Fl_Group *Elasticity3D_Menu; //GUI -> Defines the Elasticity 3D Submenu from Weak Form Menu;
  Fl_Group *ElasticTruss_Menu; //GUI -> Defines the Elastic Truss Submenu from Weak Form Menu;
  Fl_Group *PositionalTruss_Menu; //GUI -> Defines the Positional Truss Submenu from Weak Form Menu;
  Fl_Group *Stokes_Menu; //GUI -> Defines the Stokes Submenu from Weak Form Menu;
  Fl_Group *NavierStokes_Menu; //GUI -> Defines the Navier Stokes Submenu from Weak Form Menu;
  Fl_Group *Poisson_Menu; //GUI -> Defines the Poisson Submenu from Weak Form Menu;
  Fl_Group *L2Projection_Menu; //GUI -> Defines the L2 Projection Submenu from Weak Form Menu;

  Fl_Choice *PhysicalGroups; //GUI -> Defines the Physical Group menu from Weak Form Menu;
  static Fl_Menu_Item menu_PhysicalGroups[];

  Fl_Choice *WeakForm; //GUI -> Defines the Weak Form menu from Weak Form Menu;
  static Fl_Menu_Item menu_WeakForm[];

  Fl_Choice *BoundaryCondition; //GUI -> Defines the Boundary Condition menu from L2 Projection SubMenu;
  static Fl_Menu_Item menu_BoundaryCondition[];

  Fl_Value_Input *young; //GUI -> Defines the Elasticity Modulus variable intput from Weak Form Menu;
  Fl_Value_Input *poisson; //GUI -> Defines the Poisson Ratio variable intput from Weak Form Menu;
  Fl_Value_Input *area; //GUI -> Defines the Area variable intput from Weak Form Menu;
  Fl_Value_Input *density; //GUI -> Defines the Density variable intput from Weak Form Menu;
  Fl_Value_Input *viscosity; //GUI -> Defines the Viscosity variable intput from Weak Form Menu;
  Fl_Value_Input *Nsteps; //GUI -> Defines the Number of Steps variable intput from Weak Form Menu;
  Fl_Value_Input *x; //GUI -> Defines the X (val2[0]) variable intput from Weak Form Menu;
  Fl_Value_Input *y; //GUI -> Defines the Y (val2[1]) variable intput from Weak Form Menu;
  Fl_Value_Input *z; //GUI -> Defines the Z (val2[2]) variable intput from Weak Form Menu;
  Fl_Check_Button *plane_stress; //GUI -> Defines the Plane Stress variable intput from Weak Form Menu (0 for Plane Strain; 1 for Plane Stress);

  Fl_Button *Button_Apply_WeakForm; //GUI -> Defines the apply button for the WeakForm Menu widget.

  /* ========================= CREATION FUNCTIONS ========================= */

public:

  void fWF(Fl_Double_Window *Window); //GUI -> Defines the Weak Form Menu function. Also, constains all Weak Form Submenus;

  void fElasticity2D(Fl_Group *group); //GUI -> Defines the Elasticity 2D Submenu function;
  void fElasticityPositional2D(Fl_Group *group); //GUI -> Defines the Elasticity Positional 2D Submenu function;
  void fElasticity3D(Fl_Group *group); //GUI -> Defines the Elasticity 3D Submenu function;
  void fElasticTruss(Fl_Group *group); //GUI -> Defines the Elastic Truss Submenu function;
  void fPositionalTruss(Fl_Group *group); //GUI -> Defines the Positional Truss Submenu function;
  void fStokes(Fl_Group *group); //GUI -> Defines the Stokes Submenu function;
  void fNavierStokes(Fl_Group *group); //GUI -> Defines the Navier Stokes Submenu function;
  void fPoisson(Fl_Group *group); //GUI -> Defines the Poisson Submenu function;
  void fL2Projection(Fl_Group *group); //GUI -> Defines the L2 Projection Submenu function.

  /* ========================= CALLBACK FUNCTIONS ========================= */

  void fApply_WF_cb();

  void fElasticity2D_cb(); //GUI -> Defines the Elasticity 2D Callback Submenu function;
  void fElasticityPositional2D_cb(); //GUI -> Defines the Elasticity Positional 2D Callback Submenu function;
  void fElasticity3D_cb(); //GUI -> Defines the Elasticity 3D Callback Submenu function;
  void fElasticTruss_cb(); //GUI -> Defines the Elastic Truss Callback Submenu function;
  void fPositionalTruss_cb(); //GUI -> Defines the Positional Truss Callback Submenu function;
  void fStokes_cb(); //GUI -> Defines the Stokes Callback Submenu function;
  void fNavierStokes_cb(); //GUI -> Defines the Navier Stokes Callback Submenu function;
  void fPoisson_cb(); //GUI -> Defines the Poisson Callback Submenu function;
  void fL2Projection_cb(); //GUI -> Defines the L2 Projection Callback Submenu function.

private:

  inline void fInline_Apply_WF(Fl_Button*, void*);
  static void fStatic_Apply_WF(Fl_Button*, void*);

  inline void fInline_Elasticity2D(Fl_Menu_*, void*);
  static void fStatic_Elasticity2D(Fl_Menu_*, void*);

  inline void fInline_ElasticityPositional2D(Fl_Menu_*, void*);
  static void fStatic_ElasticityPositional2D(Fl_Menu_*, void*);

  inline void fInline_Elasticity3D(Fl_Menu_*, void*);
  static void fStatic_Elasticity3D(Fl_Menu_*, void*);

  inline void fInline_ElasticTruss(Fl_Menu_*, void*);
  static void fStatic_ElasticTruss(Fl_Menu_*, void*);

  inline void fInline_PositionalTruss(Fl_Menu_*, void*);
  static void fStatic_PositionalTruss(Fl_Menu_*, void*);

  inline void fInline_Stokes(Fl_Menu_*, void*);
  static void fStatic_Stokes(Fl_Menu_*, void*);

  inline void fInline_NavierStokes(Fl_Menu_*, void*);
  static void fStatic_NavierStokes(Fl_Menu_*, void*);

  inline void fInline_Poisson(Fl_Menu_*, void*);
  static void fStatic_Poisson(Fl_Menu_*, void*);

  inline void fInline_L2Projection(Fl_Menu_*, void*);
  static void fStatic_L2Projection(Fl_Menu_*, void*);


//************************************
//ANALYSIS MENU COMPONENTS - FSArl GUI
//************************************

  /* ========================= WIDGETS ========================= */

public:

  Fl_Group *Analysis_Menu; //GUI -> Defines the Analysis Menu;
  Fl_Choice *AnalysisType; //GUI -> Defines the Analysis choice from Analysis Menu;
  static Fl_Menu_Item menu_Analysis[];

  Fl_Choice *SolverType; //GUI -> Defines the Solver Type choice from Analysis Menu;
  static Fl_Menu_Item menu_SolverType[];

  Fl_Group *NonLinear_Menu; //GUI -> Defines the Non Linear Submenu form Analysis Menu;
  Fl_Group *Increm_Transient_Menu; //GUI -> Defines the Incremental and Transient Submenu form Analysis Menu;

  Fl_Choice *CompMesh; //GUI -> Defines the Comp Mesh choice from Analysis Menu;
  static Fl_Menu_Item menu_CompMesh[];

  //As Steps Number variable, the Fl_Value_Input Nsteps will be used.
  Fl_Value_Input *tolerance; //GUI -> Defines the tolerance variable intput from Analysis Menu;
  Fl_Value_Input *maxInterations; //GUI -> Defines the Maximum Interations variable intput from Analysis Menu;

  Fl_Button *Apply_Analysis; //GUI -> Defines the apply button for the Analysis Menu widget.


  /* ========================= CREATION FUNCTIONS ========================= */

public:

  void fAnalysis(Fl_Double_Window *window); //GUI -> Defines the Analysis Menu function;
  void fNonLinear(Fl_Group *group); //GUI -> Defines the NonLinear Submenu function;
  void fIncrem_Transient(Fl_Group *group); //GUI -> Defines the Incremental and Transient Submenu function. Since this two have the same variables, they have the same structure too.

  /* ========================= CALLBACK FUNCTIONS ========================= */

public:

  void fLinear_cb(); //GUI -> Defines the Linear Callback Submenu function;
  void fNonLinear_cb(); //GUI -> Defines the Non Linear Callback Submenu function;
  void fIncrem_Transient_cb(); //GUI -> Defines the Incremental and Transient Callback Submenu function;
  void fApply_Analysis_cb(); //GUI -> Defines the Apply Analysis button Callback function;

private:

  inline void fInline_Linear(Fl_Menu_*, void*);
  static void fStatic_Linear(Fl_Menu_*, void*);

  inline void fInline_NonLinear(Fl_Menu_*, void*);
  static void fStatic_NonLinear(Fl_Menu_*, void*);

  inline void fInline_Increm_Transient(Fl_Menu_*, void*);
  static void fStatic_Increm_Transient(Fl_Menu_*, void*);

  inline void fInline_Apply_Analysis(Fl_Button*, void*);
  static void fStatic_Apply_Analysis(Fl_Button*, void*);


//******************************************
//POS-PROCESSING MENU COMPONENTS - FSArl GUI
//******************************************

  /* ========================= WIDGETS ========================= */

public:

  Fl_Group *PProcessing_Menu; //GUI -> Defines the Pos-Processing Menu;
  Fl_Output *vtk; //GUI -> Defines the  vtk_path output from Pos-Processing Menu.

  /* ========================= CREATION FUNCTIONS ========================= */

public:

  void fPProcessing(Fl_Double_Window *window); //GUI -> Defines the Pos-Processing Menu function;

  /* ========================= CALLBACK FUNCTIONS ========================= */

public:

  void fInput_vtk_cb(); //GUI -> Defines the functionality from vtk Fl_Output;

private:

  inline void fInline_Input_vtk(Fl_Output*, void*);
  static void fStatic_Input_vtk(Fl_Output*, void*);


//**************************************
//TOP HEADER MENU COMPONENTS - FSArl GUI
//**************************************

  /* ========================= WIDGETS ========================= */

public:

  Fl_Menu_Bar *TopHeader_Menu; //GUI -> Defines the Top Header Menu;
  static Fl_Menu_Item Topheader_Menu[]; //GUI -> Defines the itens from the Top Header Menu.

  /* ========================= CREATION FUNCTIONS ========================= */

public:

  void fTopHeader_Menu(Fl_Double_Window *window); //GUI -> Defines the Top Header Menu function.

  /* ========================= CALLBACK FUNCTIONS ========================= */

public:

  void fOpen_File_cb(); //GUI -> Defines the Open File option callback from File item from Top Header Menu;
  void fNew_File_cb(); //GUI -> Defines the New File option callback from File item from Top Header Menu.

private:

  inline void fInline_Open_File(Fl_Menu_*, void*);
  static void fStatic_Open_File(Fl_Menu_*, void*);

  inline void fInline_New_File(Fl_Menu_*, void*);
  static void fStatic_New_File(Fl_Menu_*, void*);


//************************************
//TOOL BAR MENU COMPONENTS - FSArl GUI
//************************************

  /* ========================= WIDGETS ========================= */

public:

  Fl_Group *ToolBar_Menu; //GUI -> Defines the Tool Bar Menu.
  Fl_Button *Button_Geometry; //GUI -> Tool that shows the Geometry Menu;
  Fl_Button *Button_WeakForm; //GUI -> Tool that shows the Weak Form Menu;
  Fl_Button *Button_Analysis; //GUI -> Tool that show the Analysis Menu
  Fl_Button *Button_PosProcressing; //GUI -> Tool that show the Pos-Processing Menu

  Fl_PNG_Image *png_image_G; //GUI -> Defines the image from the Geometry button from Tool Bar Menu;
  Fl_PNG_Image *png_image_WF_act; //GUI -> Defines the image from the Weak Form activated button from Tool Bar Menu;
  Fl_PNG_Image *png_image_WF_desact; //GUI -> Defines the image from the Weak Form desactivated button from Tool Bar Menu;
  Fl_PNG_Image *png_image_A; //GUI -> Defines the image from the Analysis button from Tool Bar Menu;
  Fl_PNG_Image *png_image_PP; //GUI -> Defines the image from the Pos-Processing button from Tool Bar Menu.

  /* ========================= CREATION FUNCTIONS ========================= */

public:

  void fToolBar(Fl_Double_Window *window); //GUI -> Defines the Tool Bar creation function;
  void fButton_Geometry(Fl_Group *group); //GUI -> Defines the Geometry Button creation function;
  void fButton_WeakForm(Fl_Group *group); //GUI -> Defines the Weak Form Button creation function;
  void fButton_Analysis(Fl_Group *group); //GUI -> Defines the Analysis Button creation function;
  void fButton_PosProcressing(Fl_Group *group); //GUI -> Defines the Pos-Processing Button creation function.
  
  /* ========================= CALLBACK FUNCTIONS ========================= */

public:

  void fButton_Geometry_cb(); //GUI -> Defines the Geometry Button callback from Tool Bar Menu;
  void fButton_WeakForm_cb(); //GUI -> Defines the Weak Form Button callback from Tool Bar Menu;
  void fButton_Analysis_cb(); //GUI -> Defines the Analysis Button callback from Tool Bar Menu;
  void fButton_PosProcressing_cb(); //GUI -> Defines the Pos-Processing Button callback from Tool Bar Menu.

private:

  inline void fInline_Button_Geometry(Fl_Button*, void*);
  static void fStatic_Button_Geometry(Fl_Button*, void*);

  inline void fInline_Button_WeakForm(Fl_Button*, void*);
  static void fStatic_Button_WeakForm(Fl_Button*, void*);

  inline void fInline_Button_Analysis(Fl_Button*, void*);
  static void fStatic_Button_Analysis(Fl_Button*, void*);

  inline void fInline_Button_PProcessing(Fl_Button*, void*);
  static void fStatic_Button_PPospricessing(Fl_Button*, void*);

};
#endif
