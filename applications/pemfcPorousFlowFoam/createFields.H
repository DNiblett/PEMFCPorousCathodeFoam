Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, simple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());


singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

#include "createMRF.H"
#include "createFvOptions.H"

Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);


Info<< "Reading Permeability/Porosity Fields\n" << endl;

dimensionedScalar poreDiameter("poreDiameter", dimLength, transportProperties);
dimensionedScalar rho("rho", dimDensity, transportProperties);
dimensionedScalar nu("nu", dimViscosity, transportProperties);

volScalarField porosity
(
    IOobject
    (
        "porosity",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Calculating Permeability Field \n" << endl;
volScalarField permeability = pow(porosity, 3) * pow(poreDiameter, 2) / (150 * (1 - porosity) * (1 - porosity));

// Optional directional permeabilities allow the GDL resistance to represent
// the strong through-plane/in-plane anisotropy of carbon-paper substrates and
// the additional compaction beneath flow-field lands.  When the fields are
// absent, the Carman--Kozeny permeability above is retained in both directions
// so existing cases remain backwards compatible.
Info<< "Reading optional directional permeability fields\n" << endl;
volScalarField permeabilityInPlane
(
    IOobject
    (
        "permeabilityInPlane",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    permeability
);

volScalarField permeabilityThroughPlane
(
    IOobject
    (
        "permeabilityThroughPlane",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    permeability
);
Info<< "Calculating Forchheimer Field \n" << endl;
volScalarField C_F = 1.75 / ( sqrt(pow(porosity, 3.0) * 150));
// Compute Forchheimer coefficient field
Info<< "Min/max in-plane permeability " << min(permeabilityInPlane).value()
    << ' ' << max(permeabilityInPlane).value() << endl;
Info<< "Min/max through-plane permeability " << min(permeabilityThroughPlane).value()
    << ' ' << max(permeabilityThroughPlane).value() << endl;
