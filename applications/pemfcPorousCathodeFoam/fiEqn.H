{

/*--------------------------------------------------------------------------------
Solver - fuelCellFoam
File - fiEqn.H
Author - Dr Daniel Niblett, Newcastle University, October 2024.
Co-authors - ...

Follows inspiration from A. Colli et al.

Updates the boundary value of fi (inside the membrane/electrolyte phase) - it is not resolved in the air-space mesh. It also solves for cell potential, knowing the applied/consumed averaged cell current.

1. First solve the Tafel/Butler-Volmer Equation on the catalyst surface j_k = j_0c*exp((fi_s - fi_cathode - E0)/b_c)
2. The ionic flux from the membrane takes the anode potential in the membrane j_ohm = k/delta * (fi_cathode - fi_anode)
3. Equate the two equations for local current conservation: f(fi_cathode) = j_0c*exp((fi_s - fi_cathode - E0)/b_c) - k/delta * (fi_cathode - fi_anode)
4. Solve f(f_cathode) using Netwon-Raphson iterations: fi_cathode^(n+1) = fi_cathode^(n) - relaxation*(j_k - j_ohm)/(-j_k/b - j/delta)
5. Check mass conservation and update the shared fi_s (overall cell potential) based on the applied current density:
6. fi_s^(n+1) = fi_s^(n) - relaxation* (sum(j_k * cellArea) - appliedCurrent)/( sum(j_k * cellArea)/b_c)

Notes: V_anode = 0;
Notes: V_cell = V_cathode - V_anode;
Notes: negative kinetic properties for cathode reaction;

----------------------------------------------------------------------------------*/

// Read the applied/consumed current
scalar I = readScalar(controlProperties.lookup("I"));

scalar relaxationFactor = readScalar(controlProperties.lookup("relaxationFactor"));
word operatingMode = controlProperties.lookupOrDefault<word>
(
    "operatingMode",
    "galvanostatic"
);
const bool potentiostatic = (operatingMode == "potentiostatic");
const bool galvanostatic = (operatingMode == "galvanostatic");
if (!potentiostatic && !galvanostatic)
{
    FatalErrorInFunction
        << "Unsupported operatingMode '" << operatingMode
        << "'. Use galvanostatic or potentiostatic." << exit(FatalError);
}

// Read the fuel cell properties for membrane and kinetics
scalar C_0 = readScalar(transportProperties.lookup("C_0")); // reference concentration 
scalar b_c = readScalar(transportProperties.lookup("b_c")); // cathode tafel slope
scalar j0_c = readScalar(transportProperties.lookup("j_0c")); // anode exchange current density
scalar b_a = readScalar(transportProperties.lookup("b_a")); // anode tafel slope
scalar j0_a = readScalar(transportProperties.lookup("j_0a")); // anode exchange current density
scalar delta_m = readScalar(transportProperties.lookup("delta_m")); // anode tafel slope
scalar k_m = readScalar(transportProperties.lookup("k_m")); // anode exchange current density
scalar Vcell = readScalar(Results.lookup("Vcell")); // get the current predicted cell potential (V)
const scalar VcellSetpoint = controlProperties.lookupOrDefault<scalar>
(
    "VcellSetpoint",
    Vcell
);
const scalar potentiostaticCurrentTolerance =
    controlProperties.lookupOrDefault<scalar>
    (
        "potentiostaticCurrentTolerance",
        1e-7
    );
scalar E0_c = readScalar(transportProperties.lookup("E0_c")); // cathode tafel slope
scalar F = readScalar(transportProperties.lookup("F")); // Faraday Constant (C/mol)
scalar R = readScalar(transportProperties.lookup("R")); // ideal gas constant
scalar T = readScalar(transportProperties.lookup("T")); // temperature K
scalar alpha_a = readScalar(transportProperties.lookup("alpha_a")); // temperature K
scalar alpha_c = readScalar(transportProperties.lookup("alpha_c")); // temperature K

// Read the kinetic model chosen (Tafel or Butler Volmer)
word kineticModel = word(transportProperties.lookup("kineticModel")); 


//Vcell = E0_c - Vcell;

if (potentiostatic)
{
    // The terminal voltage is prescribed; the total current is solved below.
    Vcell = VcellSetpoint;
    I = controlProperties.lookupOrDefault<scalar>
    (
        "initialCurrent",
        max(readScalar(Results.lookup("I")), scalar(SMALL))
    );
}

       if (Vcell < 0)
       {
        Vcell = 0.1; // Enforce the minimum threshold
       }
       
       if (Vcell > 3)
       {Vcell = 1; // Enforce the max threshold
       }
      

// extract patch locations for cathode wall  
label patchi = mesh.boundaryMesh().findPatchID("catalyst"); 

//calculate the overall cell current density
scalar totalArea = 0.0;
        {
        	totalArea = gSum
                (
                	mesh.magSf().boundaryField()[patchi]
                );
        }
        
scalar currentDensity = I/totalArea;

//calculate the anode overpotential fi_a
scalar fi_as = 0; // anode solid phase potential at 0;
scalar eta_a = 0;
scalar fi_a = 0;
#include "updateAnodePotential.H"


//load oxygen concentration field on the catalyst wall boundary
const volScalarField& C_O2 = mesh.lookupObject<const volScalarField>("C_O2");
const scalarField& C_O2_face = C_O2.boundaryField()[patchi];
const scalarField& kineticMultiplierFace = kineticMultiplier.boundaryField()[patchi];
const scalarField& contactResistanceFace = contactResistance.boundaryField()[patchi];
if (min(kineticMultiplierFace) <= SMALL || min(contactResistanceFace) < -SMALL)
{
    FatalErrorInFunction
        << "kineticMultiplier must be positive and contactResistance non-negative"
        << exit(FatalError);
}
const scalarField j0cFace(j0_c*kineticMultiplierFace);
// The local resistance is evaluated explicitly from the preceding catalyst
// iterate.  This is a stable Picard linearisation of j*R_contact and is exact
// for the legacy R_contact=0 case.
scalarField Vkin(Vcell - Flux.boundaryFieldRef()[patchi]*contactResistanceFace);
forAll(Vkin, facei)
{
    Vkin[facei] = max(min(Vkin[facei], scalar(1.5)), scalar(-0.5));
}

// Calculate the local current density (Flux) using (C_O2) - (j0 and bc are negative for fuel cell)       
// Tafel
Flux.boundaryFieldRef()[patchi] = j0cFace * (C_O2_face/C_0)*exp((Vkin -fi.boundaryFieldRef()[patchi] - E0_c)/b_c); // kinetic flux

Info << "Calculate Tafel flux cathode is: " << gAverage(Flux.boundaryFieldRef()[patchi]) << endl;
// Butler-Volmer
scalar B_cf = (R*T)/(alpha_c*F);        //forwards
scalar B_cb = (R*T)/((1-alpha_c)*F);   //backwards
Flux.boundaryFieldRef()[patchi] = j0cFace * (C_O2_face/C_0) * exp((Vkin -fi.boundaryFieldRef()[patchi] - E0_c)/B_cf) - j0cFace * exp((Vkin -fi.boundaryFieldRef()[patchi] - E0_c)/B_cb); // kinetic flux
Info << "Calculate BV flux cathode is: " << gAverage(Flux.boundaryFieldRef()[patchi]) << endl;


// Initialise the Netwon-Raphson Iterative Loop
int N = 1000;
scalar tolerance = 1e-7;
//scalar minC = 1e-6;

for (int iter =0; iter < N; ++iter)
{
Vkin = Vcell - Flux.boundaryFieldRef()[patchi]*contactResistanceFace;
forAll(Vkin, facei)
{
    Vkin[facei] = max(min(Vkin[facei], scalar(1.5)), scalar(-0.5));
}
//Info << "Iteration: " << iter << endl;

//Info<< "Min/max fi:" << min(fi.boundaryField()[patchi]) << ' '
 //   << max(fi.boundaryField()[patchi]) << endl;
    
//C_O2_face = (C_O2_face/C_O2_face)*C_0;
if (kineticModel == "ButlerVolmer")
{
// START - Butler-Volmer for Cathode
const scalarField j_k_bv = j0cFace * (C_O2_face/C_0) * exp((Vkin -fi.boundaryFieldRef()[patchi] - E0_c)/B_cf) 
	                 - j0cFace * exp((Vkin -fi.boundaryFieldRef()[patchi] - E0_c)/B_cb);
	                 
const scalarField j_ohm_bv = (k_m/delta_m)*(fi.boundaryFieldRef()[patchi]- fi_a);
const scalarField f_bv = j_k_bv - j_ohm_bv;
const scalarField df_bv = j0cFace * (C_O2_face/C_0) * (-1/B_cf) *exp((Vkin -fi.boundaryFieldRef()[patchi] - E0_c)/B_cf) 
	                 - j0cFace * (-1/B_cb) * exp((Vkin -fi.boundaryFieldRef()[patchi] - E0_c)/B_cb)
	                 - (k_m/delta_m);
	                             
//fi.boundaryFieldRef()[patchi] = fi.boundaryFieldRef()[patchi] - relaxationFactor*(f_bv/df_bv);       
fi.boundaryFieldRef()[patchi] -= relaxationFactor * (f_bv / df_bv);     
	                 
//Info << "Calculate BV fi cathode is: " << gAverage(fi_test) << endl;

// Update flux Butler-Volmer
Flux.boundaryFieldRef()[patchi] = j0cFace * (C_O2_face/C_0) * exp((Vkin -fi.boundaryFieldRef()[patchi] - E0_c)/B_cf) 
	                 - j0cFace * exp((Vkin -fi.boundaryFieldRef()[patchi] - E0_c)/B_cb);
}
// END - Butler-Volmer for Cathode
	                 

// START - Tafel for Cathode
if (kineticModel == "Tafel") 
{
const scalarField j_k = j0cFace*(C_O2_face/C_0)*exp((Vkin -fi.boundaryFieldRef()[patchi] - E0_c)/b_c);
const scalarField j_ohm = (k_m/delta_m)*(fi.boundaryFieldRef()[patchi]- fi_a);
const scalarField f = j_k - j_ohm;
const scalarField df = -(j_k/b_c) - (k_m/delta_m);

// Newton Iterate
fi.boundaryFieldRef()[patchi] = fi.boundaryFieldRef()[patchi] - relaxationFactor*(f/df);

// Update Flux Tafel
Flux.boundaryFieldRef()[patchi] = -j0cFace*(C_O2_face/C_0)*exp((Vkin-fi.boundaryFieldRef()[patchi]-E0_c)/b_c);
}
// END - Tafel for Cathode

// Calculate the current on the electrode (catalyst)
scalar predictedCurrent = 0.0; // totalCurrent at cathode
        {
        	predictedCurrent = gSum
                (
                	mesh.magSf().boundaryField()[patchi]
                      * Flux.boundaryFieldRef()[patchi]
                );
        }

//Info << "predicted current is : " << predictedCurrent << endl;
//Info << "actual current is : " << I << endl;

if (galvanostatic && kineticModel == "ButlerVolmer")
{
// Newton iterate (BV)
scalar dfpredictedCurrent = 0.0; // totalCurrent at cathode
        {
        	dfpredictedCurrent = gSum
                (
                	mesh.magSf().boundaryField()[patchi]
			      * (j0cFace * (C_O2_face/C_0) * (1/B_cf) *exp((Vkin -fi.boundaryFieldRef()[patchi] - E0_c)/B_cf) 
	                 - j0cFace * (1/B_cb) * exp((Vkin -fi.boundaryFieldRef()[patchi] - E0_c)/B_cb))
                );
        }

scalar Vcell_old = Vcell;
scalar f_Vcell = -predictedCurrent - I;
scalar df_Vcell = -dfpredictedCurrent;


Vcell = Vcell - relaxationFactor*(f_Vcell/df_Vcell);
}

if (galvanostatic && kineticModel == "Tafel")
{
Vcell = Vcell - relaxationFactor*((predictedCurrent - I))/(predictedCurrent/b_c);
}
       
// Info << "predicted current is : " << -predictedCurrent << endl;
// Info << "actual current is : " << I << endl;
// Info << "Cell Potential is : " << Vcell << endl;


if (potentiostatic)
{
    // Fixed Vcell: update the anode loss with the cathode-integrated current,
    // then repeat the coupled catalyst/anode solve at the same voltage.
    const scalar Inew = max(-predictedCurrent, scalar(SMALL));
    const scalar currentResidual = mag(Inew - I);
    I = Inew;
    currentDensity = I/totalArea;
    #include "updateAnodePotential.H"

    if (iter > 0 && currentResidual < potentiostaticCurrentTolerance)
    {
        Info << "Potentiostatic current converged at " << I
            << " A in " << iter << " iterations, residual = "
            << currentResidual << endl;
        electrochemConverged = true;
        break;
    }
}
// Then keep looping around - break criteria for galvanostatic operation.
else if (mag(-predictedCurrent - I) < tolerance)
        {
        scalar residual_V = mag(-predictedCurrent - I);
        scalar eta_c = Vcell - gAverage(fi.boundaryFieldRef()[patchi]) - E0_c;
            Info << "Cathode potential converged at: " << eta_c << " V in " << iter << " iterations, " << "Final residual = " << residual_V << endl;
            electrochemConverged = true;
            break;
        }     
}

if (kineticModel == "Tafel")
{
// Update Flux (Tafel)
Vkin = Vcell - Flux.boundaryFieldRef()[patchi]*contactResistanceFace;
forAll(Vkin, facei)
{
    Vkin[facei] = max(min(Vkin[facei], scalar(1.5)), scalar(-0.5));
}
Flux.boundaryFieldRef()[patchi] = -j0cFace * (C_O2_face/C_0)*exp((Vkin -fi.boundaryFieldRef()[patchi] -E0_c)/b_c); // kinetic flux
}

if (kineticModel == "ButlerVolmer")
{
// Update Flux (BV)
Vkin = Vcell - Flux.boundaryFieldRef()[patchi]*contactResistanceFace;
forAll(Vkin, facei)
{
    Vkin[facei] = max(min(Vkin[facei], scalar(1.5)), scalar(-0.5));
}
Flux.boundaryFieldRef()[patchi] =-( j0cFace * (C_O2_face/C_0) * exp((Vkin -fi.boundaryFieldRef()[patchi] - E0_c)/B_cf) 
	                 - j0cFace * exp((Vkin -fi.boundaryFieldRef()[patchi] - E0_c)/B_cb));
}

scalar predictedCurrent = 0.0; // totalCurrent at cathode
        {
        	predictedCurrent = gSum
                (
                	mesh.magSf().boundaryField()[patchi]
                      * Flux.boundaryField()[patchi]
                );
        }

fi.boundaryFieldRef()[patchi] = fi.boundaryFieldRef()[patchi]; // kinetic flux   
//Flux.write();
//C_O2.write();
//fi.write(); 

// correct the sign to consume oxygen (Tafel)
//Flux.boundaryFieldRef()[patchi] = -j0_c * (C_O2_face/C_0)*exp((Vcell -fi.boundaryFieldRef()[patchi] -E0_c)/b_c); // kinetic flux
// Butler Volmer
//Flux.boundaryFieldRef()[patchi] =- (j0_c * (C_O2_face/C_0) * exp((Vcell -fi.boundaryFieldRef()[patchi] - E0_c)/B_cf) 
//	                 - j0_c * exp((Vcell -fi.boundaryFieldRef()[patchi] - E0_c)/B_cb));
        
// update cell potential and current density for reporting.
scalar cellOverpotential = Vcell;

 if (cellOverpotential < 0)
       {
        cellOverpotential = 0; // Enforce the minimum threshold
       }
       
scalar currentDensity2 = predictedCurrent/totalArea;

//fi.correctBoundaryConditions();

// update the values of Vcell and I in the constant folder
Results.set("Vcell", cellOverpotential);
Results.set("I", predictedCurrent);

// update eta
scalar eta_c = gAverage(Vkin - fi.boundaryFieldRef()[patchi]) - E0_c;
scalar eta_m = gAverage(fi_a - fi.boundaryFieldRef()[patchi]);

// report values in the log

Info<< "Total current predicted = " << predictedCurrent << " A "<< endl;
Info<< "Total current applied = " << I << " A "<< endl;
Info<< "Predicted current density  = " << currentDensity2 << " A/m2 "<< endl;
Info<< "Catalyst area = " << totalArea << " m2 "<< endl;
Info<< "Min/max fi:" << min(fi).value() << ' '
    << max(fi).value() << endl;
Info<< "Overpotentials: eta_anode = " << eta_a << " V, eta_cathode = "<< eta_c << " V, eta_membrane = " << eta_m << " V" << endl;
Info<< "Equalibrium Potential = " << E0_c << " V "<< endl;  
Info<< "Anode Membrane Potential = " << fi_a << " V "<< endl;
Info<< "Cathode Membrane Potential = " << gAverage(fi.boundaryField()[patchi]) << " V "<< endl;
Info<< "Cell Voltage = " << cellOverpotential << " V "<< endl;       

}

// ************************************************************************* //
