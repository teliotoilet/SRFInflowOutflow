/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

Class
    Foam::SRFInflowOutflowFvPatchVectorField

Description
    Freestream velocity condition to be used in conjunction with the single
    rotating frame (SRF) model (see: SRFModel class)

    Given the free stream velocity in the absolute frame, the condition
    applies the appropriate rotation transformation in time and space to
    determine the local velocity using:

        \f[
            U_p = cos(\theta)*U_{Inf} + sin(theta) (n^UInf) - U_{p,srf}
        \f]

    where
    \vartable
        U_p     = patch velocity [m/s]
        U_{Inf} = free stream velocity in the absolute frame [m/s]
        theta   = swept angle [rad]
        n       = axis direction of the SRF
        U_{p,srf} = SRF velocity of the patch
    \endvartable


    \heading Patch usage

    \table
        Property        | Description               | Required | Default value
        UInf            | freestream velocity       | yes      |
        relative        | UInf relative to the SRF? | no       |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            SRFInflowOutflowVelocity;
        UInf            uniform (0 0 0);
        relative        no;
        value           uniform (0 0 0);    // initial value
    }
    \endverbatim

SeeAlso
    Foam::SRFInflowOutflowFvPatchVectorField.C

SourceFiles
    SRFInflowOutflowFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef SRFInflowOutflowFvPatchVectorField_H
#define SRFInflowOutflowFvPatchVectorField_H

#include "inletOutletFvPatchFields.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class SRFInflowOutflowFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class SRFInflowOutflowFvPatchVectorField
:
    public inletOutletFvPatchVectorField
{
    // Private data

        //- Is the supplied inlet value relative to the SRF
        Switch relative_;

        //- Velocity of the free stream in the absolute frame [m/s]
        vector UInf_;


public:

    //- Runtime type information
    TypeName("SRFInflowOutflow");


    // Constructors

        //- Construct from patch and internal field
        SRFInflowOutflowFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        SRFInflowOutflowFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given SRFInflowOutflowFvPatchVectorField
        //  onto a new patch
        SRFInflowOutflowFvPatchVectorField
        (
            const SRFInflowOutflowFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        SRFInflowOutflowFvPatchVectorField
        (
            const SRFInflowOutflowFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new SRFInflowOutflowFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        SRFInflowOutflowFvPatchVectorField
        (
            const SRFInflowOutflowFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new SRFInflowOutflowFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the velocity at infinity
            const vector& UInf() const
            {
                return UInf_;
            }

            //- Return reference to the velocity at infinity to allow adjustment
            vector& UInf()
            {
                return UInf_;
            }


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
