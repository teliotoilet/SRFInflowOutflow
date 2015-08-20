/*---------------------------------------------------------------------------* \
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

\*---------------------------------------------------------------------------*/

#include "SRFInflowOutflowFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "SRFModel.H"
#include "steadyStateDdtScheme.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SRFInflowOutflowFvPatchVectorField::
SRFInflowOutflowFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    inletOutletFvPatchVectorField(p, iF),
    relative_(false),
    UInf_(vector::zero),
    calcFrac_(true)
{}


Foam::SRFInflowOutflowFvPatchVectorField::
SRFInflowOutflowFvPatchVectorField
(
    const SRFInflowOutflowFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    inletOutletFvPatchVectorField(ptf, p, iF, mapper),
    relative_(ptf.relative_),
    UInf_(ptf.UInf_),
    calcFrac_(ptf.calcFrac_)
{}


Foam::SRFInflowOutflowFvPatchVectorField::
SRFInflowOutflowFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    inletOutletFvPatchVectorField(p, iF),
    relative_(dict.lookupOrDefault("relative", false)),
    UInf_(dict.lookup("UInf")),
    calcFrac_(dict.lookupOrDefault("calcFraction", false))
{
    this->phiName_ = dict.lookupOrDefault<word>("phi","phi");

    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
}


Foam::SRFInflowOutflowFvPatchVectorField::
SRFInflowOutflowFvPatchVectorField
(
    const SRFInflowOutflowFvPatchVectorField& srfvpvf
)
:
    inletOutletFvPatchVectorField(srfvpvf),
    relative_(srfvpvf.relative_),
    UInf_(srfvpvf.UInf_),
    calcFrac_(srfvpvf.calcFrac_)
{}


Foam::SRFInflowOutflowFvPatchVectorField::
SRFInflowOutflowFvPatchVectorField
(
    const SRFInflowOutflowFvPatchVectorField& srfvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    inletOutletFvPatchVectorField(srfvpvf, iF),
    relative_(srfvpvf.relative_),
    UInf_(srfvpvf.UInf_),
    calcFrac_(srfvpvf.calcFrac_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::SRFInflowOutflowFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Get reference to the SRF model
    const SRF::SRFModel& srf =
        db().lookupObject<SRF::SRFModel>("SRFProperties");


    word ddtScheme
    (
        this->dimensionedInternalField().mesh()
       .ddtScheme(this->dimensionedInternalField().name())
    );

    if (ddtScheme == fv::steadyStateDdtScheme<scalar>::typeName)
    {
        // If not relative to the SRF include the effect of the SRF
        if (!relative_)
        {
            refValue() = UInf_ - srf.velocity(patch().Cf());

            vector omg( srf.omega().value() );       // note: r = [x,y,z]
            vector gradu(       0, -omg[2],  omg[1] ); // grad( (Omega x r) . i )
            vector gradv(  omg[2],       0, -omg[0] ); // grad( (Omega x r) . j )
            vector gradw( -omg[1],  omg[0],       0 ); // grad( (Omega x r) . k )
            forAll( *this, faceI )
            {
                refGrad()[faceI].component(0) = patch().Sf()[faceI]/patch().magSf()[faceI] & gradu;
                refGrad()[faceI].component(1) = patch().Sf()[faceI]/patch().magSf()[faceI] & gradv;
                refGrad()[faceI].component(2) = patch().Sf()[faceI]/patch().magSf()[faceI] & gradw;
            }
        }
        // If already relative to the SRF simply supply the inlet value
        // as a fixed value
        else
        {
            refValue() = UInf_;
            refGrad() = vector::zero;
        }
    }
    else
    {
//        scalar time = this->db().time().value();
//        scalar theta = time*mag(srf.omega().value());
//
//        refValue() =
//            cos(theta)*UInf_ + sin(theta)*(srf.axis() ^ UInf_)
//          - srf.velocity(patch().Cf());
        refValue() = UInf_ - srf.velocity(patch().Cf());

        vector omg( srf.omega().value() );       // note: r = [x,y,z]
        vector gradu(       0, -omg[2],  omg[1] ); // grad( (Omega x r) . i )
        vector gradv(  omg[2],       0, -omg[0] ); // grad( (Omega x r) . j )
        vector gradw( -omg[1],  omg[0],       0 ); // grad( (Omega x r) . k )
        // this gives uniform 0...
        //refGrad().component(0) = patch().Sf() & gradu;
        //refGrad().component(1) = patch().Sf() & gradv;
        //refGrad().component(2) = patch().Sf() & gradw;
        forAll( *this, faceI )
        {
            refGrad()[faceI].component(0) = patch().Sf()[faceI]/patch().magSf()[faceI] & gradu;
            refGrad()[faceI].component(1) = patch().Sf()[faceI]/patch().magSf()[faceI] & gradv;
            refGrad()[faceI].component(2) = patch().Sf()[faceI]/patch().magSf()[faceI] & gradw;

            // note: this makes the BC select inflow/outflow based on initial
            //       rather than instantanteous values
            valueFraction()[faceI] = 1.0 - pos(UInf_ & patch().Sf()[faceI]);
        }

    }

    // Set the inlet-outlet choice based on the direction of the freestream
    // note: this causes a switching between inflow/outflow on the outer boundary...
    if( calcFrac_ )
    {
        valueFraction() = 1.0 - pos(refValue() & patch().Sf());
    }
    else valueFraction() = 0.0;

    ////////////////////////////////////////////////////////////
    // DEBUG ///////////////////////////////////////////////////
    // Notes:
    // * Sf points outward from computational domain
//    forAll(*this, faceI)
//    {
//        vector bf(patch().Cf()[faceI]);
//        Info<< patch().name() << " " 
//            << bf << ":"
//            << " value" << refValue()[faceI]
//            << " grad" << refGrad()[faceI]
//            << " dx=" << 1/patch().deltaCoeffs()[faceI]
//            << " frac=" << valueFraction()[faceI] << endl;
//    }
    ////////////////////////////////////////////////////////////

    mixedFvPatchField<vector>::updateCoeffs();
}


void Foam::SRFInflowOutflowFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("relative") << relative_ << token::END_STATEMENT << nl;
    os.writeKeyword("UInf") << UInf_ << token::END_STATEMENT << nl;
    os.writeKeyword("calculateFraction") << calcFrac_ << token::END_STATEMENT << nl;
    os.writeKeyword("phi") << this->phiName_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        SRFInflowOutflowFvPatchVectorField
    );
}


// ************************************************************************* //
