
file(GLOB integration 
    IntegrationQuadrature.h
    IntegrationQuadrature.cpp
    IntegrationQuadrature11.h
    IntegrationQuadrature11.cpp
    DomainIntegration.h
    DomainIntegration.cpp
    IntRule.h
    IntRule.cpp
    IntRule0d.h
    IntRule0d.cpp
    IntRule1d.h
    IntRule1d.cpp
    IntRuleHexahedron.h
    IntRuleHexahedron.cpp
    IntRuleQuad.h
    IntRuleQuad.cpp
    IntRuleTetrahedron.h
    IntRuleTetrahedron.cpp
    IntRuleTriangle.h
    IntRuleTriangle.cpp
)

set(INTEGRAL ${integration} PARENT_SCOPE)


