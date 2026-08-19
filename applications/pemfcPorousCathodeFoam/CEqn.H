 {
        // concentration fields
	dictionary speciesSolveDict = mesh.solutionDict().subDict("solvers").subDict("C_species");
	// The mask is zero in fixed-grid topology ``solid'' pixels.  Applying it
	// to both terms is critical: a high Darcy penalty alone makes velocity
	// small, but it does not mathematically prevent concentration leakage.
	volScalarField speciesAccess(porosity*speciesTransportMask);
	// Linear interpolation alone would give a non-zero face coefficient between
	// a 1-mask and 0-mask cell.  This boolean face gate is one only when both
	// neighbouring cells are accessible, making a topology solid an exact
	// no-advection/no-diffusion barrier rather than merely a high resistance.
	surfaceScalarField speciesFaceAccess
	(
	    pos(fvc::interpolate(speciesTransportMask) - scalar(0.999))
	);
	surfaceScalarField phi_vol_avg = phi * fvc::interpolate(porosity)
	    * speciesFaceAccess;

        forAll(speciesNameList, i)
        {
        
        // V2: face-normal projection of a diagonal effective-diffusivity
        // tensor.  This preserves the legacy model when both tortuosities
        // match, while permitting independent in-plane and through-plane O2
        // transport without modifying the mesh or boundary conditions.
        volScalarField D_inPlane
        (
            speciesList[i].D()*speciesAccess*oxygenDiffusivityMultiplier/tauInPlane
        );
        volScalarField D_throughPlane
        (
            speciesList[i].D()*speciesAccess
           *oxygenDiffusivityThroughPlaneMultiplier/tauThroughPlane
        );
        surfaceScalarField nZ
        (
            mesh.Sf().component(vector::Z)/mesh.magSf()
        );
        surfaceScalarField D_eff
        (
            speciesFaceAccess*
            (
                fvc::interpolate(D_inPlane)*(scalar(1) - sqr(nZ))
              + fvc::interpolate(D_throughPlane)*sqr(nZ)
            )
        );
		dimensionedScalar solidCellAnchor
		(
		    "solidCellAnchor", dimless/dimTime, 1e6
		);
        

            fvScalarMatrix CEqn
            (
                //  fvm::ddt(speciesList[i])
                fvm::div(phi_vol_avg, speciesList[i], "div(phi,C_species)") // remove convection for now
                - fvm::laplacian(D_eff, speciesList[i], "laplacian(D,C_species)") //  speciesList[i].D()+Dturb      
			+ fvm::Sp((scalar(1) - speciesTransportMask)*solidCellAnchor, speciesList[i])
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
