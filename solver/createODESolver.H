IOdictionary odeSolver
(
    IOobject
    (
        "odeSolver",
        runTime.system(),
        runTime,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

autoPtr<ODESolver> solver = ODESolver::New(pointKineticsEquations, odeSolver);
