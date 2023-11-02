
file(GLOB integration 
    BoundaryIntegrationQuadrature.h
    BoundaryIntegrationQuadrature.cpp
    IntegrationQuadrature.h
    IntegrationQuadrature.cpp
    IntegrationQuadrature11.h
    IntegrationQuadrature11.cpp
    DomainIntegration.h
    DomainIntegration.cpp
)

set(INTEGRAL ${integration} PARENT_SCOPE)


