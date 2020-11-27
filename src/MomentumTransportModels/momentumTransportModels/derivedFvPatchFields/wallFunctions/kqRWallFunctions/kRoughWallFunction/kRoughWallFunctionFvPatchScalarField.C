/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2020 OpenFOAM Foundation
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

#include "kRoughWallFunctionFvPatchScalarField.H"
#include "nutWallFunctionFvPatchScalarField.H"
#include "momentumTransportModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

kRoughWallFunctionFvPatchScalarField::kRoughWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF),
    Ceps2_(1.9),
    Ks_(p.size(), 0.0)
{}


kRoughWallFunctionFvPatchScalarField::kRoughWallFunctionFvPatchScalarField
(
    const kRoughWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper),
    Ceps2_(ptf.Ceps2_),
    Ks_(ptf.Ks_)
{}


kRoughWallFunctionFvPatchScalarField::kRoughWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict),
    Ceps2_(dict.lookupOrDefault<scalar>("Ceps2", 1.9)),
    Ks_("Ks", dict, p.size())
{}


kRoughWallFunctionFvPatchScalarField::kRoughWallFunctionFvPatchScalarField
(
    const kRoughWallFunctionFvPatchScalarField& kwfpsf
)
:
    fixedValueFvPatchField<scalar>(kwfpsf),
    Ceps2_(kwfpsf.Ceps2_),
    Ks_(kwfpsf.Ks_)
{}


kRoughWallFunctionFvPatchScalarField::kRoughWallFunctionFvPatchScalarField
(
    const kRoughWallFunctionFvPatchScalarField& kwfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(kwfpsf, iF),
    Ceps2_(kwfpsf.Ceps2_),
    Ks_(kwfpsf.Ks_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void kRoughWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const label patchi = patch().index();

    const momentumTransportModel& turbModel =
        db().lookupObject<momentumTransportModel>
        (
            IOobject::groupName
            (
                momentumTransportModel::typeName,
                internalField().group()
            )
        );

    const nutWallFunctionFvPatchScalarField& nutw =
        nutWallFunctionFvPatchScalarField::nutw(turbModel, patchi);

    const scalarField& y = turbModel.y()[patchi];

    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();

    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    const scalar Cmu25 = pow025(nutw.Cmu());

    const scalar beta5 = pow(0.075, 0.5);
    
    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    const scalarField magGradUw(mag(Uw.snGrad())); 
		
    scalarField& kw = *this;

    // Set k wall values
    forAll(kw, facei)
    {
        //label celli = patch().faceCells()[facei];

        //scalar uTau = Cmu25*sqrt(k[celli]);

        //scalar yPlus = uTau*y[facei]/nuw[facei];

	/*
        if (yPlus > nutw.yPlusLam())
        {
            scalar Ck = -0.416;
            scalar Bk = 8.366;
            kw[facei] = Ck/nutw.kappa()*log(yPlus) + Bk;
        }
        else
        {
            scalar C = 11.0;
            scalar Cf = (1.0/sqr(yPlus + C) + 2.0*yPlus/pow3(C) - 1.0/sqr(C));
            kw[facei] = 2400.0/sqr(Ceps2_)*Cf;
        }

        kw[facei] *= sqr(uTau);
    	*/
    	
    	//Knopp Roughness modification
    	scalar ut = sqrt( (nuw[facei]+nutw[facei]) * magGradUw[facei] );
	scalar KsPlus = Ks_[facei]*ut / nuw[facei];
	kw[facei] = min(1.,KsPlus/90.)*ut*ut/beta5;
    }

    // Limit kw to avoid failure of the turbulence model due to division by kw
    kw = max(kw, small);

    fixedValueFvPatchField<scalar>::updateCoeffs();

    // TODO: perform averaging for cells sharing more than one boundary face
}


void kRoughWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    writeEntry(os, "Ceps2", Ceps2_);
    writeEntry(os, "Ks", Ks_);
    fixedValueFvPatchField<scalar>::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    kRoughWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
