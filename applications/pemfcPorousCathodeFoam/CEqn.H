 {
        // concentration fields
	dictionary speciesSolveDict = mesh.solutionDict().subDict("solvers").subDict("C_species");

        forAll(speciesNameList, i)
        {
        
        // U is the intrinsic pore velocity in the species equation, hence
        // epsilon U is the superficial molar-convective velocity.  This is
        // the conservative steady volume-averaged form in the GDL and
        // reduces exactly to the free-flow equation where epsilon = 1.
        surfaceScalarField phiPorous(phi*fvc::interpolate(porosity));

        // The resistance of a discontinuous channel--GDL interface is
        // represented with a harmonic normal-tensor diffusivity.  This keeps
        // the total diffusive face flux conservative and avoids arithmetic
        // averaging across the large bulk/GDL diffusivity jump.
        surfaceVectorField nHat(mesh.Sf()/mesh.magSf());
        surfaceSymmTensorField inverseDface(fvc::interpolate(inv(D_O2)));
        surfaceScalarField D_eff(scalar(1)/(nHat & (inverseDface & nHat)));
        

            fvScalarMatrix CEqn
            (
                //  fvm::ddt(speciesList[i])
                fvm::div(phiPorous, speciesList[i], "div(phiPorous,C_species)")
                - fvm::laplacian(D_eff, speciesList[i], "laplacian(D_O2,C_species)")
		==
		  fvOptions(speciesList[i])
	    );



            CEqn.relax();
	    fvOptions.constrain(CEqn);
	    
	    CEqn.solve(speciesSolveDict);
	    fvOptions.correct(speciesList[i]);

            // A concentration cannot be negative.  Enforce the physical
            // lower bound after the segregated solve so a locally depleted
            // cell cannot destabilise the nonlinear catalyst update.
            speciesList[i].max
            (
                dimensionedScalar
                (
                    "Cmin",
                    speciesList[i].dimensions(),
                    SMALL
                )
            );

// output to log min and max concentration of species
	    Info<< "Min/max C_" << speciesNameList[i] << ":" << min(speciesList[i]).value() << ' '
            << max(speciesList[i]).value() << endl;
	}
}

 
 
 
 
 // Previous code
 /*
 {
 fvScalarMatrix CEqn
            (
               // fvm::ddt(T)
               //fvm::div(phi, C)
              - fvm::laplacian(DC, C)
           
            );

            CEqn.relax();
            fvOptions.constrain(CEqn);
            CEqn.solve();
            fvOptions.correct(C);
            
            }
            
            */
