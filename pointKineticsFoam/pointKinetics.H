/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
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

Class
    Foam::pointKinetics

Description

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef pointKinetics_H
#define pointKinetics_H

#include "ODESystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class pointKinetics Declaration
\*---------------------------------------------------------------------------*/

class pointKinetics
:
    public ODESystem
{
    // Private data

        //- Description of data_
        scalar reactivity_;//$
        scalar totalDelayedPrecursorFraction_;
        scalar firstDelayedPrecursorFraction_;
        scalar firstDelayedPrecursorDecayConstant_;
        scalar neutronGenerationTime_;//s

    // Private Member Functions

        //- Disallow default bitwise copy construct
        pointKinetics(const pointKinetics&);

        //- Disallow default bitwise assignment
        void operator=(const pointKinetics&);


public:

    // Constructors

        //- Construct from components
        pointKinetics(
            const scalar reactivity,
            const scalar totalDelayedPrecursorFraction,
//          const scalar firstDelayedPrecursorFraction,
            const scalar firstDelayedPrecursorDecayConstant,
            const scalar neutronGenerationTime
        );

    // Member Functions

        //- Return the number of equations in the system
        label nEqns() const;

        //- Calculate the drivatives in dydx
        void derivatives
        (
            const scalar x,
            const scalarField& y,
            scalarField& dydx
        ) const;

        //- Calculate the Jacobian of the system
        //  Need by the stiff-system solver
        void jacobian
        (
            const scalar x,
            const scalarField& y,
            scalarField& dfdx,
            scalarSquareMatrix& dfdy
        ) const;

        //- Reactivity feedback
        void reactivityFeedback(scalar reactivity);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
