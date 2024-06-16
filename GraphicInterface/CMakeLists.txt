
file(GLOB graphicinterface 
    WindowConstructor.h
    OpenGLWindow/OpenGL_WC.h
    OpenGLWindow/OpenGL_WC.cpp
    OpenFLTKWindow/OpenFLTKWindow.cpp
    WeakFormMenu/Elasticity2DMenu.cpp
    WeakFormMenu/ElasticityPositional2DMenu.cpp
    WeakFormMenu/Elasticity3DMenu.cpp
    WeakFormMenu/ElasticTrussMenu.cpp
    WeakFormMenu/PositionalTrussMenu.cpp
    WeakFormMenu/StokesDMenu.cpp
    WeakFormMenu/NavierStokesMenu.cpp
    WeakFormMenu/PoissonMenu.cpp
    WeakFormMenu/L2Projection2DMenu.cpp
    WeakFormMenu/WeakFormMenu.cpp
    WeakFormMenu/GeometryMenu.cpp
    WeakFormMenu/OpenFoamMenu.cpp
    TopHeaderMenu/TopHeaderMenu.cpp
    TopHeaderMenu/TopHeaderFunctions.cpp
    ToolBarMenu/Button_Analysis.cpp
    ToolBarMenu/Button_Geometry.cpp
    ToolBarMenu/Button_WeakForm.cpp
    ToolBarMenu/Button_PosProcessing.cpp
    ToolBarMenu/ToolBarMenu.cpp
    GeometryMenu/GeometryFunctions.cpp
    GeometryMenu/Geometry_Inl_Sta_Functions.cpp
    WeakFormMenu/StokesMenu.cpp
    WeakFormMenu/WeakFormFunctions.cpp
    WeakFormMenu/WeakForm_Inl_Sta_Functions.cpp
    ToolBarMenu/ToolBar_Inl_Sta_Functions.cpp
    ToolBarMenu/ToolBarFunctions.cpp
    ScriptDisplay/ScriptPopup.cpp
    ScriptDisplay/ScriptDisplay.cpp
    ScriptDisplay/Script_Inl_Sta_Functions.cpp
    ScriptDisplay/ScriptFunctions.cpp
    TopHeaderMenu/TopHeader_Inl_Sta_Functions.cpp
    AnalysisMenu/Analysis_Inl_Sta_Functions.cpp
    AnalysisMenu/AnalysisFunctions.cpp
    AnalysisMenu/AnalysisMenu.cpp
    AnalysisMenu/NonLinearMenu.cpp
    AnalysisMenu/Increment&TransientMenu.cpp
    PosProcessingMenu/PProcessingMenu.cpp
    PosProcessingMenu/PProcessingFunction.cpp
    PosProcessingMenu/PProcessing_Inl_Sta_Functions.cpp
    
    WindowConstructor.cpp
)

set(GRAPHICINTERFACE ${graphicinterface} PARENT_SCOPE)


