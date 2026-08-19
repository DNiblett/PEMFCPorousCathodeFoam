    // Report Permeability Field
    //#include "transportProperties.H" // Include the file that has your transport properties

// Read transport properties

//word useForchheimer = word(transportProperties.lookup("useForchheimer"));
//scalar C_F = readScalar(transportProperties.lookup("C_F")); // Faraday Constant (C/mol)

// If packedBed

   
    //+ fvc::div((C_F / sqrt(permeability)) * mag(U) * U)
    // Momentum predictor

    MRF.correctBoundaryVelocity(U);

    // Darcy resistance is evaluated from the full symmetric permeability
    // tensor.  Its isotropic trace is implicit and the anisotropic remainder
    // is explicit, so the converged equation is nu*inv(K).U without assuming
    // that the principal directions align with x/y/z.
    volSymmTensorField inversePermeability(inv(permeability));
    volScalarField darcyImplicit(nu*tr(inversePermeability)/scalar(3));
    volVectorField darcyTensorCorrection
    (
        IOobject
        (
            "darcyTensorCorrection",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector
        (
            "zero",
            U.dimensions()/dimTime,
            vector::zero
        )
    );
    darcyTensorCorrection =
        nu*(inversePermeability & U) - darcyImplicit*U;

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::div(phi, U)
      + MRF.DDt(U)
      + turbulence->divDevReff(U)
      + fvm::Sp(forchheimerCoefficient*mag(U), U)
      + fvm::Sp(darcyImplicit, U)
      + darcyTensorCorrection
     ==
        fvOptions(U)
    );
    
     

    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (simple.momentumPredictor())
    {
        solve(UEqn == -fvc::grad(p));

        fvOptions.correct(U);
    }
