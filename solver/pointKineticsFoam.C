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
#include "ODESolver.H"
#include "IOdictionary.H"
#include "Time.H"

#include "pointKinetics.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList args(argc, argv);

    #include "createTime.H"
    #include "createPointKineticsEquations.H"
    #include "createScalarFields.H"
    #include "createODESolver.H"

    while (runTime.running())
    {
        scalar t = runTime.timeOutputValue();
        scalar deltaT = runTime.deltaT().value();
        
        pointKineticsEquations.derivatives(t, y, dy);
        solver->solve(t, y, deltaT);
        
        runTime.setDeltaT(deltaT);
        runTime.setTime(t, runTime.timeIndex() + 1);
    }

    return 0;
}
// ************************************************************************* //
