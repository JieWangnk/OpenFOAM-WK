/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
     \\/     M anipulation  |
---------------------------------------------------------------------------*/

#include "modularWKPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

namespace Foam
{

// Constructors

modularWKPressureFvPatchScalarField::modularWKPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    order_(readLabel(dict.lookup("order"))),
    R_(readScalar(dict.lookup("R"))),
    C_(readScalar(dict.lookup("C"))),
    Z_(readScalar(dict.lookup("Z"))),
    p1_(0.0), // Will be calculated in updateCoeffs
    p0_(readScalar(dict.lookup("p0"))),
    p_1_(dict.lookupOrDefault("p_1", p0_)),
    q0_(0.0), // Will be calculated in updateCoeffs
    q_1_(readScalar(dict.lookup("q_1"))),
    q_2_(dict.lookupOrDefault("q_2", q_1_)),
    q_3_(dict.lookupOrDefault("q_3", q_2_))
{
    // Set the initial pressure value of the patch from p0
    fixedValueFvPatchScalarField::operator==(p0_);
}


modularWKPressureFvPatchScalarField::modularWKPressureFvPatchScalarField
(
    const modularWKPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    order_(ptf.order_),
    R_(ptf.R_),
    C_(ptf.C_),
    Z_(ptf.Z_),
    p1_(ptf.p1_),
    p0_(ptf.p0_),
    p_1_(ptf.p_1_),
    q0_(ptf.q0_),
    q_1_(ptf.q_1_),
    q_2_(ptf.q_2_),
    q_3_(ptf.q_3_)
{}

modularWKPressureFvPatchScalarField::modularWKPressureFvPatchScalarField
(
    const modularWKPressureFvPatchScalarField& fvmpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(fvmpsf, iF),
    phiName_(fvmpsf.phiName_),
    order_(fvmpsf.order_),
    R_(fvmpsf.R_),
    C_(fvmpsf.C_),
    Z_(fvmpsf.Z_),
    p1_(fvmpsf.p1_),
    p0_(fvmpsf.p0_),
    p_1_(fvmpsf.p_1_),
    q0_(fvmpsf.q0_),
    q_1_(fvmpsf.q_1_),
    q_2_(fvmpsf.q_2_),
    q_3_(fvmpsf.q_3_)
{}


// Member Functions

void modularWKPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // --- 1. Get the flux from the previous timestep's result ---
    const surfaceScalarField& phi =
        db().lookupObject<surfaceScalarField>(phiName_);

    // Sum the flux over the patch to get the flow rate Q
    q0_ = sum(phi.boundaryField()[this->patch().index()]);

    // --- 2. Solve the Windkessel ODE for the new pressure ---
    const scalar dt = db().time().deltaTValue();
    scalar Q_source = 0.0;
    scalar Pgrad_part = 0.0;
    scalar Pdenom = 1.0;

    switch (order_)
    {
        case 1:
            Q_source = (q0_/C_)*(1 + Z_/R_) + (Z_/dt)*(q0_ - q_1_);
            Pgrad_part = -p0_/dt;
            Pdenom = 1/dt + 1/(R_*C_);
            break;

        case 2:
            Q_source = (q0_/C_)*(1 + Z_/R_) + (Z_/dt)*(1.5*q0_ - 2*q_1_ + 0.5*q_2_);
            Pgrad_part = (-2*p0_ + 0.5*p_1_)/dt;
            Pdenom = 1.5/dt + 1/(R_*C_);
            break;

        case 3:
            Q_source = (q0_/C_)*(1 + Z_/R_) + (Z_/dt)*((11.0/6.0)*q0_ - 3*q_1_ + 1.5*q_2_ - (1.0/3.0)*q_3_);
            Pgrad_part = (-3*p0_ + 1.5*p_1_ - (1.0/3.0)*p_1_)/dt;
            Pdenom = (11.0/6.0)/dt + 1/(R_*C_);
            break;

        default:
            // Default to 1st order if an invalid order is given
            Q_source = (q0_/C_)*(1 + Z_/R_) + (Z_/dt)*(q0_ - q_1_);
            Pgrad_part = -p0_/dt;
            Pdenom = 1/dt + 1/(R_*C_);
            break;
    }

    p1_ = (Q_source - Pgrad_part) / Pdenom;


    // --- 3. Set the boundary condition value for this timestep ---
    this->operator==(p1_);


    // --- 4. Update historical values for the next timestep ---
    q_3_ = q_2_;
    q_2_ = q_1_;
    q_1_ = q0_;

    p_1_ = p0_;
    p0_ = p1_;

    fixedValueFvPatchScalarField::updateCoeffs();
}


void modularWKPressureFvPatchScalarField::write(Ostream& os) const
{
    // Use the base class to write the "type" and "value" entries
    fixedValueFvPatchScalarField::write(os);

    // Write the parameters
    os.writeKeyword("phi") << phiName_ << token::END_STATEMENT << nl;
    os.writeKeyword("order") << order_ << token::END_STATEMENT << nl;
    os.writeKeyword("R") << R_ << token::END_STATEMENT << nl;
    os.writeKeyword("C") << C_ << token::END_STATEMENT << nl;
    os.writeKeyword("Z") << Z_ << token::END_STATEMENT << nl;

    // Write the historical state for robust restarts
    os.writeKeyword("p0") << p0_ << token::END_STATEMENT << nl;
    os.writeKeyword("p_1") << p_1_ << token::END_STATEMENT << nl;
    os.writeKeyword("q_1") << q_1_ << token::END_STATEMENT << nl;
    os.writeKeyword("q_2") << q_2_ << token::END_STATEMENT << nl;
    os.writeKeyword("q_3") << q_3_ << token::END_STATEMENT << nl;
}

} // End namespace Foam

#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    makePatchTypeField(fvPatchScalarField, modularWKPressureFvPatchScalarField);
}