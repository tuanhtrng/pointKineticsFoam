/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Description

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IOmanip.H"
#include "IOdictionary.H"
#include "pointKinetics.H"
#include "ODESolver.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::validArgs.append("ODESolver");
    argList args(argc, argv);

    scalar net_reactivity = 0.2*0.0065;
    const scalar delayed_neutron_fraction = 0.0065;
    const scalar prompt_neutron_generation_time = 0.0001;//s
    const scalar delayed_neutron_decay_constant = 0.07741;//s^-1

    // initial conditions
    const scalar rho_0 = net_reactivity;//$
    const scalar y0_0 = 10.0;//MW
    const scalar y1_0 = 
        (delayed_neutron_fraction - rho_0)
       /prompt_neutron_generation_time
       /delayed_neutron_decay_constant
       *y0_0;//MW 
    const label n = 1000;
    const scalar endTime = 100.0;

    // Create the ODE system
    pointKinetics ode
    (
        net_reactivity,
        delayed_neutron_fraction,
        delayed_neutron_decay_constant,
        prompt_neutron_generation_time
    );

    dictionary dict;
    dict.add("solver", args[1]);

    // Create the selected ODE system solver
    autoPtr<ODESolver> odeSolver = ODESolver::New(ode, dict);

    // Initialise the ODE system fields

    // start time
    scalar xStart = 0.0;

    // time step
    const scalar dx = endTime/n;
    // initial displacement and velocity
    scalarField yStart(ode.nEqns());
    yStart[0] = y0_0;
    yStart[1] = y1_0;
    // integration initial step
    scalar dxEst = 0.1;
    scalar xEnd = 0.0;

    // required to store dydx
    scalarField dyStart(ode.nEqns());

    Info<< "Time"
        << "\t" << "Neutron Density (MW)"
        << endl;

    // integration loop
    for (label i=0; i<n; i++)
    {
        Info<< xStart 
            << "\t" << yStart[0]
            << endl;
        xEnd = xStart + dx;
        ode.derivatives(xStart, yStart, dyStart);
        odeSolver->solve(xStart, xEnd, yStart, dxEst);
        xStart = xEnd;
    }

    Info<< "\nEND\n" << endl;

    return 0;
}
// ************************************************************************* //
