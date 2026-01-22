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
    // Two-parameter control with backward compatibility
    // If betaT specified, use it; otherwise fall back to legacy beta
    betaT_
    (
        dict.found("betaT")
            ? dict.lookup<scalar>("betaT")
            : dict.lookupOrDefault<scalar>("beta", 0.2)
    ),
    betaN_(dict.lookupOrDefault<scalar>("betaN", 0.0)),
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
    betaT_(ptf.betaT_),
    betaN_(ptf.betaN_),
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
    betaT_(swvf.betaT_),
    betaN_(swvf.betaN_),
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
        // Two-parameter directional stabilization using directionMixed approach
        //
        // Formula: valueFraction = backflowMask * (betaN*sqr(n) + betaT*(I-sqr(n)))
        //
        // Where:
        //   sqr(n) = n⊗n = normal projection tensor
        //   (I - sqr(n)) = tangential projection tensor
        //
        // Physical interpretation:
        //   - betaT controls tangential backflow suppression (vortices)
        //   - betaN controls normal backflow suppression (flow reversal)
        //   - For Windkessel compatibility, use betaN=0 (normal responds to pressure)
        //
        // This is matrix-coupled at each PIMPLE outer iteration for proper p-U coupling

        // Clamp effective betas to [0,1] for safety
        const scalar effBetaT = min(max(betaT_ * dampingFactor_, 0.0), 1.0);
        const scalar effBetaN = min(max(betaN_ * dampingFactor_, 0.0), 1.0);

        // Backflow mask: 1 if phi<0 (backflow), 0 otherwise
        const scalarField backflowMask(pos0(-phip));

        // Normal vector field and projection tensors
        const vectorField n = patch().nf();
        const symmTensorField normalProj = sqr(n);        // n⊗n
        const symmTensorField tangProj = I - normalProj;  // I - n⊗n

        // Combined valueFraction for two-parameter control
        valueFraction() = backflowMask * (effBetaN * normalProj + effBetaT * tangProj);
    }

    // refValue stays at zero (target for backflow suppression)
    // refGrad stays at zero (zeroGradient behavior)

    // Call base class - NO evaluate() call after this!
    // Let solver framework call evaluate() at the right time in PIMPLE loop
    directionMixedFvPatchVectorField::updateCoeffs();
}


void Foam::stabilizedWindkesselVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);

    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);

    // Write two-parameter control
    writeEntry(os, "betaT", betaT_);
    writeEntry(os, "betaN", betaN_);

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
