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

    Inherits from directionMixedFvPatchVectorField to provide:
    - Proper implicit matrix treatment for PIMPLE algorithm convergence
    - Directional control: only tangential velocity is damped during backflow
    - Normal velocity responds freely to pressure (Windkessel compatibility)

    This approach is physically motivated for pressure-driven boundaries:
    - The pressure BC (Windkessel) should drive the normal velocity
    - Only tangential backflow (vortices, recirculation) should be suppressed

    Behavior:
    - Outflow (phi >= 0): zeroGradient for all velocity components
    - Inflow/backflow (phi < 0):
        * Normal component: zeroGradient (responds to pressure)
        * Tangential components: damped toward zero

    Usage:
    \verbatim
    outlet
    {
        type                  stabilizedWindkesselVelocity;
        phi                   phi;
        beta                  0.5;      // Damping strength (0-1)
        enableStabilization   true;
        dampingFactor         1.0;      // Additional damping control
        value                 uniform (0 0 0);
    }
    \endverbatim

    Parameters:
    - phi: Name of flux field (default: "phi")
    - beta: Stabilization coefficient (0-1, default: 0.5)
      Higher = more tangential backflow suppression
    - enableStabilization: Toggle stabilization (default: true)
    - dampingFactor: Additional damping multiplier (default: 1.0)

See also
    Foam::directionMixedFvPatchVectorField
    Foam::pressureInletOutletVelocityFvPatchVectorField

\*---------------------------------------------------------------------------*/

#ifndef stabilizedWindkesselVelocityFvPatchVectorField_H
#define stabilizedWindkesselVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "directionMixedFvPatchFields.H"

namespace Foam
{

class stabilizedWindkesselVelocityFvPatchVectorField
:
    public directionMixedFvPatchVectorField
{
    // Private Data

        //- Name of flux field
        word phiName_;

        //- Stabilization coefficient (0-1)
        //  Controls how much tangential backflow is suppressed
        //  0 = no stabilization (pure zeroGradient)
        //  1 = full suppression (tangential velocity -> 0 for backflow)
        scalar beta_;

        //- Flag to enable/disable stabilization
        bool enableStabilization_;

        //- Additional damping factor (0-1)
        //  Effective damping = beta * dampingFactor
        scalar dampingFactor_;

        //- Fluid density (for future extensions)
        scalar rho_;


public:

    //- Runtime type information
    TypeName("stabilizedWindkesselVelocity");


    // Constructors

        //- Construct from patch, internal field and dictionary
        stabilizedWindkesselVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given field onto a new patch
        stabilizedWindkesselVelocityFvPatchVectorField
        (
            const stabilizedWindkesselVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fieldMapper&
        );

        //- Construct as copy setting internal field reference
        stabilizedWindkesselVelocityFvPatchVectorField
        (
            const stabilizedWindkesselVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<vector>> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<vector>>
            (
                new stabilizedWindkesselVelocityFvPatchVectorField(*this, iF)
            );
        }


    //- Destructor
    virtual ~stabilizedWindkesselVelocityFvPatchVectorField() = default;


    // Member Functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


} // End namespace Foam

#endif
