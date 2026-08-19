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


Info<< "Reading porous-flow transport tensors\n" << endl;

dimensionedScalar nu("nu", dimViscosity, transportProperties);
dimensionedScalar poreDiameter("poreDiameter", dimLength, transportProperties);
const scalar permeabilityInPlaneMultiplier =
    transportProperties.lookupOrDefault<scalar>("permeabilityInPlaneMultiplier", 1.0);
const scalar permeabilityThroughPlaneMultiplier =
    transportProperties.lookupOrDefault<scalar>("permeabilityThroughPlaneMultiplier", 1.0);
dimensionedScalar forchheimerCoefficient
(
    "forchheimerCoefficient",
    dimless/dimLength,
    transportProperties
);

Info<< "Reading porosity field\n" << endl;
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

if (min(porosity).value() <= SMALL || max(porosity).value() > scalar(1) + SMALL)
{
    FatalErrorInFunction
        << "porosity must be in (0, 1]" << exit(FatalError);
}

// The unified channel/GDL momentum model approaches the free-flow limit as
// epsilon -> 1 (very large K); the GDL tensor decreases as epsilon falls.
volScalarField solidFraction(max(scalar(1) - porosity, scalar(1e-6)));
volScalarField permeabilityIsotropic
(
    pow(porosity, 3)*sqr(poreDiameter)/(scalar(150)*sqr(solidFraction))
);
volSymmTensorField permeability
(
    IOobject
    (
        "permeability",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedSymmTensor("zero", dimArea, symmTensor::zero)
);
permeability.replace(symmTensor::XX, permeabilityInPlaneMultiplier*permeabilityIsotropic);
permeability.replace(symmTensor::YY, permeabilityInPlaneMultiplier*permeabilityIsotropic);
permeability.replace(symmTensor::ZZ, permeabilityThroughPlaneMultiplier*permeabilityIsotropic);

Info<< "Permeability tensor calculated from porosity; Forchheimer coefficient = "
    << forchheimerCoefficient.value() << " m^-1" << nl << endl;
