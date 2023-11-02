
file(GLOB analysis 
    Analysis.h
    Analysis.cpp
    LinearAnalysis.h
    LinearAnalysis.cpp
    IncrementalAnalysis.h
    IncrementalAnalysis.cpp
    NonLinearAnalysis.h
    NonLinearAnalysis.cpp
    TransientAnalysis.h
    TransientAnalysis.cpp
    Assemble.h
    Assemble.cpp
)

set(ANALYSIS ${analysis} PARENT_SCOPE)


