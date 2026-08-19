    // Report Permeability Field
    //#include "transportProperties.H" // Include the file that has your transport properties

// Read transport properties

//word useForchheimer = word(transportProperties.lookup("useForchheimer"));
//scalar C_F = readScalar(transportProperties.lookup("C_F")); // Faraday Constant (C/mol)

// If packedBed

   
    //+ fvc::div((C_F / sqrt(permeability)) * mag(U) * U)
    // Momentum predictor

    MRF.correctBoundaryVelocity(U);

    // The vector matrix accepts a scalar implicit sink.  Treat the dominant
    // in-plane Darcy resistance implicitly and add the z-direction difference
    // explicitly.  At convergence this is exactly
    // nu*diag(1/Kip, 1/Kip, 1/Ktp).U, while retaining a well-conditioned
    // momentum matrix for the very small through-plane GDL velocities.
    volVectorField throughPlaneDarcyCorrection
    (
        IOobject
        (
            "throughPlaneDarcyCorrection",
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
    throughPlaneDarcyCorrection.replace
    (
        vector::Z,
        nu*(scalar(1)/permeabilityThroughPlane
          - scalar(1)/permeabilityInPlane)*U.component(vector::Z)
    );

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::div(phi, U)
      + MRF.DDt(U)
      + turbulence->divDevReff(U)
      // The Forchheimer contribution is scalar; use K_parallel because the
      // channel-aligned velocity dominates its magnitude in this continuum
      // model.  The directional Darcy part below supplies K_perp exactly at
      // convergence.
      + fvm::Sp( (((1-porosity)*C_F) / sqrt(permeabilityInPlane)) * mag(U), U)
      + fvm::Sp(nu / permeabilityInPlane, U)
      + throughPlaneDarcyCorrection
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
