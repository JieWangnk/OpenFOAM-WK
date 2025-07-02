/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Class
    Foam::modularWKPressureFvPatchScalarField

Description
    A modular, 3-element Windkessel pressure boundary condition.
    This condition is applied to the pressure field 'p' and calculates the
    outlet pressure based on the flux 'phi' from the previous timestep.
    This provides a numerically stable, lagged coupling mechanism.

\*---------------------------------------------------------------------------*/

#ifndef modularWKPressureFvPatchScalarField_H
#define modularWKPressureFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

namespace Foam
{

class modularWKPressureFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private Data

        //- Name of the flux field
        word phiName_;

        //- Finite difference order for dQ/dt
        label order_;

        //- Windkessel parameters
        scalar R_;
        scalar C_;
        scalar Z_;

        //- Historical pressure values (p1_ is current, p0_ is t-dt, etc.)
        scalar p1_; // Pressure for the current time step (calculated value)
        scalar p0_; // Pressure from previous time step
        scalar p_1_; // Pressure from t-2*dt

        //- Historical flow values (q0_ is from previous step, q_1_ is t-2*dt)
        scalar q0_;
        scalar q_1_;
        scalar q_2_;
        scalar q_3_;


public:

    //- Runtime type information
    TypeName("modularWKPressure");


    // Constructors

        //- Construct from patch and internal field and dictionary
        modularWKPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping
        modularWKPressureFvPatchScalarField
        (
            const modularWKPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        modularWKPressureFvPatchScalarField
        (
            const modularWKPressureFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

    //- Destructor
    virtual ~modularWKPressureFvPatchScalarField() = default;


    // Member Functions

        //- Update the coefficients
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};

} // End namespace Foam
#endif