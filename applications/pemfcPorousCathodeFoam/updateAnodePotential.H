// Update the anode/membrane-side potential from the current operating point.
// `currentDensity`, `fi_a` and `eta_a` are supplied by fiEqn.H.
eta_a = 0.0;
fi_a = 0.0;

if (kineticModel == "Tafel")
{
    fi_a = fi_as - b_a*std::log(currentDensity/j0_a);
}
else if (kineticModel == "ButlerVolmer")
{
    const scalar B_af = (R*T)/(alpha_a*F);
    const scalar B_ab = (R*T)/(-(1.0-alpha_a)*F);
    scalar eta_a_new = eta_a;

    for (int anodeIter = 0; anodeIter < 10000; ++anodeIter)
    {
        const scalar j_k_anode =
            j0_a*(std::exp(eta_a/B_af) - std::exp(eta_a/B_ab));
        const scalar residual = j_k_anode - currentDensity;
        const scalar derivative =
            j0_a*
            ((1.0/B_af)*std::exp(eta_a/B_af)
           - (1.0/B_ab)*std::exp(eta_a/B_ab));
        eta_a_new = eta_a - relaxationFactor*(residual/derivative);

        if (mag(eta_a_new - eta_a) < 1e-7)
        {
            break;
        }
        eta_a = eta_a_new;
    }
    fi_a = fi_as - eta_a;
}
else
{
    FatalErrorInFunction
        << "Unsupported kineticModel '" << kineticModel
        << "'. Use Tafel or ButlerVolmer." << exit(FatalError);
}
