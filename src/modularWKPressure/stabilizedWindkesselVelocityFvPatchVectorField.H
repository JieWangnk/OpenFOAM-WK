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
        * Normal component: controlled by betaN (0=free, 1=clamped)
        * Tangential components: controlled by betaT (0=free, 1=clamped)

    Usage:
    \verbatim
    outlet
    {
        type                  stabilizedWindkesselVelocity;
        phi                   phi;
        betaT                 0.2;      // Tangential damping (0-1)
        betaN                 0.0;      // Normal damping (0-1)
        enableStabilization   true;
        dampingFactor         1.0;      // Additional damping multiplier
        value                 uniform (0 0 0);
    }
    \endverbatim

    Parameters:
    - phi: Name of flux field (default: "phi")
    - betaT: Tangential stabilization coefficient (0-1, default: 0.2)
      Controls damping of tangential backflow (vortices, recirculation)
    - betaN: Normal stabilization coefficient (0-1, default: 0.0)
      Controls damping of normal backflow (use 0 for Windkessel compatibility)
    - beta: Legacy parameter - if specified alone, sets betaT=beta, betaN=0
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

        //- Tangential stabilization coefficient (0-1)
        //  Controls how much tangential backflow is suppressed
        //  0 = no tangential stabilization (pure zeroGradient)
        //  1 = full suppression (tangential velocity -> 0 for backflow)
        scalar betaT_;

        //- Normal stabilization coefficient (0-1)
        //  Controls how much normal backflow is suppressed
        //  0 = no normal stabilization (Windkessel-compatible, default)
        //  1 = full suppression (normal velocity -> 0 for backflow)
        //  Note: betaN > 0 may interfere with Windkessel pressure coupling
        scalar betaN_;

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
