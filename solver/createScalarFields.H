IOdictionary initialConditions
(
    IOobject
    (
        "initialConditions",
        runTime.constant(),
        runTime,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

scalarField y(pointKineticsEquations.nEqns());

y[0] = readScalar(initialConditions.lookup("neutronDensity"));
y[1] = readScalar(initialConditions.lookup("firstDelayedPrecursorDensity"));

scalarField dy(pointKineticsEquations.nEqns());