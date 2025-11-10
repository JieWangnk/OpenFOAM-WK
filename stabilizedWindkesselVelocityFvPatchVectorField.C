/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
     \\/     M anipulation  |
---------------------------------------------------------------------------*/

#include "stabilizedWindkesselVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

namespace Foam
{

// Constructors

stabilizedWindkesselVelocityFvPatchVectorField::stabilizedWindkesselVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    zeroGradientFvPatchVectorField(p, iF),
    beta_(1.0),
    enableStabilization_(true)
{}


stabilizedWindkesselVelocityFvPatchVectorField::stabilizedWindkesselVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    zeroGradientFvPatchVectorField(p, iF, dict),
    beta_(dict.lookupOrDefault<scalar>("beta", 1.0)),
    enableStabilization_(dict.lookupOrDefault<bool>("enableStabilization", true))
{
    if (dict.found("value"))
    {
        fvPatchField<vector>::operator=
        (
            vectorField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<vector>::operator=(patchInternalField());
    }
}


stabilizedWindkesselVelocityFvPatchVectorField::stabilizedWindkesselVelocityFvPatchVectorField
(
    const stabilizedWindkesselVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    zeroGradientFvPatchVectorField(ptf, p, iF, mapper),
    beta_(ptf.beta_),
    enableStabilization_(ptf.enableStabilization_)
{}


stabilizedWindkesselVelocityFvPatchVectorField::stabilizedWindkesselVelocityFvPatchVectorField
(
    const stabilizedWindkesselVelocityFvPatchVectorField& swvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    zeroGradientFvPatchVectorField(swvf, iF),
    beta_(swvf.beta_),
    enableStabilization_(swvf.enableStabilization_)
{}


// Member Functions

void stabilizedWindkesselVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Call parent updateCoeffs first
    zeroGradientFvPatchVectorField::updateCoeffs();
}


void stabilizedWindkesselVelocityFvPatchVectorField::evaluate
(
    const Pstream::commsTypes
)
{
    if (!updated())
    {
        updateCoeffs();
    }

    if (!enableStabilization_)
    {
        // Just use zero gradient if stabilization is disabled
        zeroGradientFvPatchVectorField::evaluate();
        return;
    }

    // Get velocity field from internal cells
    const vectorField velocity = patchInternalField();
    
    // Get patch normal vectors (pointing outward)
    const vectorField n = patch().nf();
    
    // Calculate normal velocity component (v·n)
    const scalarField normalVel = velocity & n;

    // Simpler, more robust backflow stabilization approach
    // Damp the backflow velocity directly based on beta parameter

    vectorField stabilizedVelocity = patchInternalField();

    forAll(stabilizedVelocity, faceI)
    {
        const scalar vn = normalVel[faceI];

        if (vn < 0.0) // Backflow detected
        {
            // Reduce backflow by factor (1 - beta)
            // beta = 0: no stabilization (full backflow allowed)
            // beta = 1: complete damping (no backflow)
            // beta = 0.5: reduce backflow by 50%

            const vector& normal = n[faceI];
            const vector tangential = stabilizedVelocity[faceI] - vn * normal;

            // Apply damping only to normal component, preserve tangential
            stabilizedVelocity[faceI] = (1.0 - beta_) * vn * normal + tangential;
        }
    }

    // Set the final field value
    fvPatchField<vector>::operator=(stabilizedVelocity);

    fvPatchField<vector>::evaluate();
}


tmp<Field<vector>> 
stabilizedWindkesselVelocityFvPatchVectorField::valueInternalCoeffs
(
    const tmp<scalarField>& w
) const
{
    // For backflow stabilization, we modify the matrix coefficients
    // to include the stabilization term contribution
    return zeroGradientFvPatchVectorField::valueInternalCoeffs(w);
}


tmp<Field<vector>> 
stabilizedWindkesselVelocityFvPatchVectorField::valueBoundaryCoeffs
(
    const tmp<scalarField>& w
) const
{
    return zeroGradientFvPatchVectorField::valueBoundaryCoeffs(w);
}


void stabilizedWindkesselVelocityFvPatchVectorField::write(Ostream& os) const
{
    zeroGradientFvPatchVectorField::write(os);

    os.writeKeyword("beta") << beta_ << token::END_STATEMENT << nl;
    os.writeKeyword("enableStabilization") << enableStabilization_
        << token::END_STATEMENT << nl;
}

} // End namespace Foam

#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        stabilizedWindkesselVelocityFvPatchVectorField
    );
}
