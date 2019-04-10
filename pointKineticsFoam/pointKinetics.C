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

\*---------------------------------------------------------------------------*/

#include "pointKinetics.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointKinetics::pointKinetics(
    const Foam::scalar reactivity,
    const Foam::scalar totalDelayedPrecursorFraction,
//  const Foam::scalar firstDelayedPrecursorFraction,
    const Foam::scalar firstDelayedPrecursorDecayConstant,
    const Foam::scalar neutronGenerationTime
)
:
    ODESystem(),
    reactivity_(reactivity),
    totalDelayedPrecursorFraction_(totalDelayedPrecursorFraction),
//  firstDelayedPrecursorFraction_(firstDelayedPrecursorFraction),
    firstDelayedPrecursorFraction_(totalDelayedPrecursorFraction),
    firstDelayedPrecursorDecayConstant_(firstDelayedPrecursorDecayConstant),
    neutronGenerationTime_(neutronGenerationTime)
{};

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

//- Return the number of equations in the system
Foam::label
Foam::pointKinetics::nEqns() const { return 2; };

//- Calculate the drivatives in dydx
void
Foam::pointKinetics::derivatives
(
    const Foam::scalar x,
    const Foam::scalarField& y,
    Foam::scalarField& dydx
) const
{
    // Neutron density
    dydx[0] =
        // Prompt neutron generation
        (reactivity_ - totalDelayedPrecursorFraction_)
       /neutronGenerationTime_
       *y[0]
      + // Delayed precursor group one
        firstDelayedPrecursorDecayConstant_
       *y[1];

    // Delayed precursor density group one
    dydx[1] = 
        // Prompt neutron
        firstDelayedPrecursorFraction_
       /neutronGenerationTime_
       *y[0]
      - // Delayed precursor group one
        firstDelayedPrecursorDecayConstant_
       *y[1];
};

//- Calculate the Jacobian of the system
//  Need by the stiff-system solver
void 
Foam::pointKinetics::jacobian
(
    const Foam::scalar x,
    const Foam::scalarField& y,
    Foam::scalarField& dfdx,
    Foam::scalarSquareMatrix& dfdy
) const
{
    dfdx[0] = 0.0;
    dfdx[1] = 0.0;

    // Prompt neutron
    dfdy[0][0] =
        (reactivity_ - totalDelayedPrecursorFraction_)
       /neutronGenerationTime_;
    // Delayed precursor
    dfdy[0][1] = firstDelayedPrecursorDecayConstant_;

    // Delayed precursor
    dfdy[1][0] = 
        firstDelayedPrecursorFraction_/neutronGenerationTime_;
    dfdy[1][1] = - firstDelayedPrecursorDecayConstant_;
};

//- Reactivity feedback
void
Foam::pointKinetics::reactivityFeedback
(
    Foam::scalar reactivity
)
{ 
    reactivity_ = reactivity; 
};
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
