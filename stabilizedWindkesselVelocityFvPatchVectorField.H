/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Class
    Foam::stabilizedWindkesselVelocityFvPatchVectorField

Description
    A backflow-stabilized velocity boundary condition for windkessel outlets.

    When backflow occurs (v·n < 0), this BC damps the backflow velocity
    using a beta factor:
        V_stabilized = (1 - beta) * V_backflow + V_tangential

    Where:
    - beta = 0.0: No stabilization (full backflow allowed)
    - beta = 0.5: Reduces backflow by 50%
    - beta = 0.9: Reduces backflow by 90% (recommended)
    - beta = 1.0: Completely blocks backflow (smooth, not hard constraint)

    This prevents numerical instability from backflow at outlet boundaries.

\*---------------------------------------------------------------------------*/

#ifndef stabilizedWindkesselVelocityFvPatchVectorField_H
#define stabilizedWindkesselVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "zeroGradientFvPatchFields.H"

namespace Foam
{

class stabilizedWindkesselVelocityFvPatchVectorField
:
    public zeroGradientFvPatchVectorField
{
    // Private Data

        //- Stabilization coefficient (typically 0.5-1.0)
        scalar beta_;

        //- Flag to enable/disable stabilization
        bool enableStabilization_;

public:

    //- Runtime type information
    TypeName("stabilizedWindkesselVelocity");


    // Constructors

        //- Construct from patch and internal field
        stabilizedWindkesselVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch and internal field and dictionary
        stabilizedWindkesselVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping
        stabilizedWindkesselVelocityFvPatchVectorField
        (
            const stabilizedWindkesselVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        stabilizedWindkesselVelocityFvPatchVectorField
        (
            const stabilizedWindkesselVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

    //- Destructor
    virtual ~stabilizedWindkesselVelocityFvPatchVectorField() = default;


    // Member Functions

        //- Update the coefficients
        virtual void updateCoeffs();
        
        //- Evaluate the patch field
        virtual void evaluate
        (
            const Pstream::commsTypes commsType = Pstream::commsTypes::blocking
        );

        //- Return the matrix diagonal coefficients corresponding to the
        //  evaluation of the value of this patchField with given weights
        virtual tmp<Field<vector>> valueInternalCoeffs
        (
            const tmp<scalarField>&
        ) const;

        //- Return the matrix source coefficients corresponding to the
        //  evaluation of the value of this patchField with given weights
        virtual tmp<Field<vector>> valueBoundaryCoeffs
        (
            const tmp<scalarField>&
        ) const;

        //- Write
        virtual void write(Ostream&) const;
};

} // End namespace Foam

#endif
