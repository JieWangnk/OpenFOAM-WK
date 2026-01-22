/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "stabilizedWindkesselVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::stabilizedWindkesselVelocityFvPatchVectorField::
stabilizedWindkesselVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    directionMixedFvPatchVectorField(p, iF),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    beta_(dict.lookupOrDefault<scalar>("beta", 0.5)),
    enableStabilization_(dict.lookupOrDefault<bool>("enableStabilization", true)),
    dampingFactor_(dict.lookupOrDefault<scalar>("dampingFactor", 1.0)),
    rho_(dict.lookupOrDefault<scalar>("rho", 1060.0))
{
    // Set initial field value
    fvPatchVectorField::operator=
    (
        vectorField("value", iF.dimensions(), dict, p.size())
    );

    // Initialize directionMixed BC parameters
    // refValue: target velocity for tangential backflow (zero to suppress)
    refValue() = Zero;

    // refGrad: gradient for zeroGradient behavior (zero)
    refGrad() = Zero;

    // valueFraction: tensor field controlling directional behavior
    // Will be updated in updateCoeffs() based on phi
    valueFraction() = Zero;
}


Foam::stabilizedWindkesselVelocityFvPatchVectorField::
stabilizedWindkesselVelocityFvPatchVectorField
(
    const stabilizedWindkesselVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fieldMapper& mapper
)
:
    directionMixedFvPatchVectorField(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    beta_(ptf.beta_),
    enableStabilization_(ptf.enableStabilization_),
    dampingFactor_(ptf.dampingFactor_),
    rho_(ptf.rho_)
{}


Foam::stabilizedWindkesselVelocityFvPatchVectorField::
stabilizedWindkesselVelocityFvPatchVectorField
(
    const stabilizedWindkesselVelocityFvPatchVectorField& swvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedFvPatchVectorField(swvf, iF),
    phiName_(swvf.phiName_),
    beta_(swvf.beta_),
    enableStabilization_(swvf.enableStabilization_),
    dampingFactor_(swvf.dampingFactor_),
    rho_(swvf.rho_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::stabilizedWindkesselVelocityFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Get flux field
    const fvsPatchField<scalar>& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);

    if (!enableStabilization_)
    {
        // No stabilization: pure zeroGradient behavior
        // valueFraction = 0 means use gradient (which is zero)
        valueFraction() = Zero;
    }
    else
    {
        // Stabilization enabled using directionMixed approach
        //
        // The valueFraction tensor controls which velocity components are
        // constrained to refValue vs allowed to float (zeroGradient)
        //
        // Key tensor: (I - n⊗n) = tangential projection
        // - Projects velocity onto the plane tangent to the patch
        // - Normal component (n direction) is NOT constrained
        // - Tangential components ARE constrained to refValue (zero)
        //
        // This is physically correct for pressure-driven boundaries:
        // - Normal velocity should respond to pressure (Windkessel)
        // - Tangential backflow (vortices) should be suppressed
        //
        // With beta and dampingFactor for adjustable stabilization strength:
        // - effectiveDamping = 0: pure zeroGradient
        // - effectiveDamping = 1: full tangential suppression

        const scalar effectiveDamping = beta_ * dampingFactor_;

        // neg(phip) returns 1 for negative values (backflow), 0 otherwise
        // sqr(patch().nf()) is the normal projection tensor n⊗n
        // (I - sqr(n)) is the tangential projection tensor
        valueFraction() = effectiveDamping * neg(phip) * (I - sqr(patch().nf()));
    }

    // refValue stays at zero (target for tangential backflow suppression)
    // refGrad stays at zero (zeroGradient behavior)

    directionMixedFvPatchVectorField::updateCoeffs();
    directionMixedFvPatchVectorField::evaluate();
}


void Foam::stabilizedWindkesselVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);

    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    writeEntry(os, "beta", beta_);
    writeEntry(os, "enableStabilization", enableStabilization_);
    writeEntry(os, "dampingFactor", dampingFactor_);
    writeEntry(os, "rho", rho_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        stabilizedWindkesselVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
