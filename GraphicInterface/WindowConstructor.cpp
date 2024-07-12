#include "WindowConstructor.h"
#include "OpenGLWindow/OpenGL_WC.h"
#include <FL/Fl_File_Chooser.H>
#include <iostream>
#include "GeometryMenu/GeometryMenu.cpp"
#include "GeometryMenu/GeometryFunctions.cpp"
#include "GeometryMenu/Geometry_Inl_Sta_Functions.cpp"
#include "WeakFormMenu/L2ProjectionMenu.cpp"
#include "WeakFormMenu/StokesMenu.cpp"
#include "WeakFormMenu/WeakFormMenu.cpp"
#include "WeakFormMenu/WeakFormFunctions.cpp"
#include "WeakFormMenu/WeakForm_Inl_Sta_Functions.cpp"
#include "WeakFormMenu/OpenFoamMenu.cpp"
#include "ToolBarMenu/ToolBarMenu.cpp"
#include "ToolBarMenu/Button_WeakForm.cpp"
#include "ToolBarMenu/Button_Analysis.cpp"
#include "ToolBarMenu/Button_Geometry.cpp"
#include "ToolBarMenu/Button_PosProcessing.cpp"
#include "ToolBarMenu/ToolBar_Inl_Sta_Functions.cpp"
#include "ToolBarMenu/ToolBarFunctions.cpp"
#include "ScriptDisplay/ScriptPopup.cpp"
#include "ScriptDisplay/ScriptDisplay.cpp"
#include "ScriptDisplay/Script_Inl_Sta_Functions.cpp"
#include "ScriptDisplay/ScriptFunctions.cpp"
#include "TopHeaderMenu/TopHeaderMenu.cpp"
#include "TopHeaderMenu/TopHeaderFunctions.cpp"
#include "TopHeaderMenu/TopHeader_Inl_Sta_Functions.cpp"
#include "TopHeaderMenu/GenerationFunctions/HeaderFile.cpp"
#include "TopHeaderMenu/GenerationFunctions/FoamFile.cpp"
#include "TopHeaderMenu/GenerationFunctions/OFfolderGenerator.cpp"
#include "TopHeaderMenu/GenerationFunctions/pFile.cpp"
#include "TopHeaderMenu/GenerationFunctions/UFile.cpp"
#include "TopHeaderMenu/GenerationFunctions/physicalPropertiesFile.cpp"
#include "TopHeaderMenu/GenerationFunctions/momentumTransportFile.cpp"
#include "TopHeaderMenu/GenerationFunctions/RunFoam.cpp"
#include "TopHeaderMenu/GenerationFunctions/controlDictFile.cpp"
#include "AnalysisMenu/Analysis_Inl_Sta_Functions.cpp"
#include "AnalysisMenu/AnalysisFunctions.cpp"
#include "AnalysisMenu/AnalysisMenu.cpp"
#include "AnalysisMenu/NonLinearMenu.cpp"
#include "AnalysisMenu/Increment&TransientMenu.cpp"
#include "AnalysisMenu/OpenFoamAnalysisMenu.cpp"
#include "PosProcessingMenu/PProcessingMenu.cpp"
#include "PosProcessingMenu/PProcessingFunction.cpp"
#include "PosProcessingMenu/PProcessing_Inl_Sta_Functions.cpp"
#include "WeakFormMenu/OpenFoamSubMenu/PressureSubMenu.cpp"
#include "WeakFormMenu/OpenFoamSubMenu/VelocitySubMenu.cpp"


/**
 Constructor
*/
WindowConstructor::WindowConstructor() {
  //Creates the FSArl Window.
   Window = new Fl_Double_Window(1536, 793, "FSArl");
    Window->color(FL_LIGHT1);
    Window->user_data((void*)(this));
    Window->align(Fl_Align(65));
    Window->resizable(Window);

    //Just a box that divides the FSArl window for the input and output variables.
    { Fl_Group* o = new Fl_Group(1159, 34, 376, 750);
      o->box(FL_BORDER_BOX);
      o->end();
    }

    //Adds the widgets:
    fGeometry(Window);
    fWF(Window);
    fAnalysis(Window);
    fPProcessing(Window);
    fScript_Menu(Window);
    fScript_Popup_Menu(Window);
    fToolBar(Window);
    fTopHeader_Menu(Window);

    //Creates and adds the OpenGL window in to the FSArl Window.
    playback = new Playback(445, 80, 714, 703, "Playback");
    Window->add(playback);
    playback->resizable();
    playback->take_focus();
    playback->end();
    
    Window->end();
    
  
}

