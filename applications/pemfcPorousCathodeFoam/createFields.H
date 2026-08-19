// * * * * * * * * * * * * * Properties fields * * * * * * * * * * * * * * * //

Info<< "Reading speciesProperties\n" << endl;

IOdictionary speciesProperties
(
    IOobject
    (
        "speciesProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

const wordList speciesNameList(speciesProperties.toc());

PtrList<species> speciesList(speciesNameList.size());

forAll(speciesNameList, i)
{
    const word& speciesName = speciesNameList[i];
    const dictionary& speciesDict = speciesProperties.subDict(speciesName);

    speciesList.set
    (
        i,
        new species
        (
            speciesDict,
            mesh
        )
    );

}


// read velocity field
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


volScalarField Flux
    (
        IOobject
        (
            "Flux",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
	dimensionedScalar("0", dimensionSet(0,-2,0,0,0,1,0), 0.0)
	 );
	 
volScalarField fi
(
    IOobject
    (
        "fi",
        runTime.timeName(),
        mesh,
        // Reading a supplied potential field provides a physically meaningful
        // initial guess for stiff potentiostatic low-voltage calculations.
        // The default field below is retained when no restart field exists.
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);
 
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

Info<< "Reading controlProperties\n" << endl;

IOdictionary controlProperties
(
    IOobject
    (
        "controlProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

IOdictionary Results
(
    IOobject
    (
        "Results",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::AUTO_WRITE
    )
);

// added parts from fuelCellFoam

dimensionedScalar tau
(
    transportProperties.lookup("tau")
);

// V2: oxygen transport may be anisotropic.  The legacy scalar tau remains a
// backward-compatible default for both directions.  x/y are in-plane and z
// is through-plane for the supplied PEMFC geometries.
const scalar tauInPlane = transportProperties.lookupOrDefault<scalar>
(
    "tauInPlane",
    tau.value()
);
const scalar tauThroughPlane = transportProperties.lookupOrDefault<scalar>
(
    "tauThroughPlane",
    tau.value()
);
if
(
    tau.value() < scalar(1)
 || tauInPlane < scalar(1)
 || tauThroughPlane < scalar(1)
)
{
    FatalErrorInFunction
        << "tau, tauInPlane and tauThroughPlane are tortuosities and must be >= 1"
        << exit(FatalError);
}
Info<< "Oxygen tortuosity: in-plane=" << tauInPlane
    << ", through-plane=" << tauThroughPlane << nl << endl;

Info<< "Reading porosity field \n" << endl;

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

// ``speciesTransportMask`` is a strict geometric accessibility field used by
// fixed-grid topology optimisation.  A value of one retains the conventional
// porous-medium transport; zero removes both species advection and diffusion
// from that cell without requiring a singular zero permeability in the flow
// solver.  It is intentionally independent of porosity: topology walls can
// use a small finite porosity/permeability for numerically conditioned flow
// penalisation while remaining exactly impermeable to species transport.
volScalarField speciesTransportMask
(
    IOobject
    (
        "speciesTransportMask", runTime.timeName(), mesh,
        IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("one", dimless, 1.0)
);
if (min(speciesTransportMask).value() < -SMALL
 || max(speciesTransportMask).value() > scalar(1) + SMALL)
{
    FatalErrorInFunction
        << "speciesTransportMask must be bounded in [0,1]"
        << exit(FatalError);
}

// Optional V2 calibration fields.  They retain legacy behaviour when absent.
// oxygenDiffusivityMultiplier and oxygenDiffusivityThroughPlaneMultiplier
// apply to the directional effective diffusivities; kineticMultiplier scales
// local ECSA/exchange-current density; and
// contactResistance is an effective local series resistance in ohm m2 (its
// field is deliberately dimensionless because the existing kinetic model uses
// dimensioned values as scalars at the catalyst boundary).
volScalarField oxygenDiffusivityMultiplier
(
    IOobject
    (
        "oxygenDiffusivityMultiplier", runTime.timeName(), mesh,
        IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("one", dimless, 1.0)
);
volScalarField oxygenDiffusivityThroughPlaneMultiplier
(
    IOobject
    (
        "oxygenDiffusivityThroughPlaneMultiplier", runTime.timeName(), mesh,
        IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("one", dimless, 1.0)
);
volScalarField kineticMultiplier
(
    IOobject
    (
        "kineticMultiplier", runTime.timeName(), mesh,
        IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("one", dimless, 1.0)
);
volScalarField contactResistance
(
    IOobject
    (
        "contactResistance", runTime.timeName(), mesh,
        IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);

#include "createMRF.H"
#include "createFvOptions.H"
#include "createPhi.H" //uncomment if convection not included
