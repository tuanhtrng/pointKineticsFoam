IOdictionary reactorProperties
(
    IOobject
    (
        "reactorProperties",
        runTime.constant(),
        runTime,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

pointKinetics pointKineticsEquations(reactorProperties);