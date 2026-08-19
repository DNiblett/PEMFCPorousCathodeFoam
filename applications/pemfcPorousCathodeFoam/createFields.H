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

// The channel has epsilon=1 and recovers bulk O2 diffusion.  In the GDL,
// porosity-dependent, directional Bruggeman corrections generate the full
// effective diffusivity tensor from properties rather than user-supplied
// directional scalar fields.
dimensionedScalar D_O2Bulk("D_O2Bulk", dimArea/dimTime, transportProperties);
const scalar bruggemanExponentInPlane = transportProperties.lookupOrDefault<scalar>
(
    "bruggemanExponentInPlane",
    1.5
);
const scalar bruggemanExponentThroughPlane = transportProperties.lookupOrDefault<scalar>
(
    "bruggemanExponentThroughPlane",
    1.5
);
if (bruggemanExponentInPlane < 0 || bruggemanExponentThroughPlane < 0)
{
    FatalErrorInFunction
        << "Bruggeman exponents must be non-negative" << exit(FatalError);
}
volScalarField D_O2InPlane(D_O2Bulk*pow(porosity, bruggemanExponentInPlane));
volScalarField D_O2ThroughPlane(D_O2Bulk*pow(porosity, bruggemanExponentThroughPlane));
Info<< "Calculating porosity-dependent O2 diffusivity tensor\n" << endl;
volSymmTensorField D_O2
(
    IOobject
    (
        "D_O2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedSymmTensor("zero", dimArea/dimTime, symmTensor::zero)
);
D_O2.replace(symmTensor::XX, D_O2InPlane);
D_O2.replace(symmTensor::YY, D_O2InPlane);
D_O2.replace(symmTensor::ZZ, D_O2ThroughPlane);

const scalar kineticMultiplier = transportProperties.lookupOrDefault<scalar>
(
    "kineticMultiplier",
    1.0
);
const scalar contactResistance = transportProperties.lookupOrDefault<scalar>
(
    "contactResistance",
    0.0
);
if (kineticMultiplier <= SMALL || contactResistance < -SMALL)
{
    FatalErrorInFunction
        << "kineticMultiplier must be positive and contactResistance non-negative"
        << exit(FatalError);
}

#include "createMRF.H"
#include "createFvOptions.H"
#include "createPhi.H" //uncomment if convection not included
