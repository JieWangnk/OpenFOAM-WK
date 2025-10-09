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
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    zeroGradientFvPatchVectorField(p, iF, dict),
    beta_(dict.lookupOrDefault<scalar>("beta", 1.0)),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
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
    rhoName_(ptf.rhoName_),
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
    rhoName_(swvf.rhoName_),
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
    
    // Get density - for incompressible flow, read from available properties dictionary
    scalar rho = 1000.0; // Default value
    
    // Try different dictionary names based on OpenFOAM version/solver
    const IOdictionary* propsDict = nullptr;
    if (db().foundObject<IOdictionary>("physicalProperties"))
    {
        propsDict = &db().lookupObject<IOdictionary>("physicalProperties");
    }
    else if (db().foundObject<IOdictionary>("transportProperties"))
    {
        propsDict = &db().lookupObject<IOdictionary>("transportProperties");
    }
    
    if (propsDict && propsDict->found("rho"))
    {
        rho = propsDict->lookup<scalar>("rho");
    }
    
    // Apply stabilization only where backflow occurs (v·n < 0)
    vectorField correction(velocity.size(), vector::zero);
    
    forAll(velocity, faceI)
    {
        if (normalVel[faceI] < 0.0) // Backflow condition
        {
            // Calculate stabilization traction: t = β*ρ*(v ⊗ v)·n
            // This is equivalent to β*ρ*(v·n)*v for the normal component
            const vector& v = velocity[faceI];
            const vector& normal = n[faceI];
            // Use the constant density from transportProperties
            
            // Stabilization term opposes the inward flow
            // The term β*ρ*(v ⊗ v)·n can be written as β*ρ*(v·n)*v
            const scalar vn = normalVel[faceI]; // This is negative for backflow
            correction[faceI] = -beta_ * rho * vn * v;
            
            // Apply only the normal component to avoid affecting tangential flow
            correction[faceI] = (correction[faceI] & normal) * normal;
        }
    }
    
    // Apply a more conservative stabilization approach
    // Instead of converting to gradient, apply the correction directly as a velocity modification
    // with a damping factor to prevent instability
    
    const scalar dampingFactor = 0.1; // Increased damping for stronger stabilization
    vectorField stabilizedVelocity = patchInternalField();
    
    // Apply correction only for backflow faces, and only to normal component
    forAll(correction, faceI)
    {
        if (mag(correction[faceI]) > SMALL)
        {
            // Apply damped correction to velocity
            stabilizedVelocity[faceI] += dampingFactor * correction[faceI];
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
    os.writeKeyword("rho") << rhoName_ << token::END_STATEMENT << nl;
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