/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "omegaRoughWallFunctionFvPatchScalarField.H"
#include "nutWallFunctionFvPatchScalarField.H"
#include "momentumTransportModel.H"
#include "fvMatrix.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

scalar omegaRoughWallFunctionFvPatchScalarField::tolerance_ = 1e-5;

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void omegaRoughWallFunctionFvPatchScalarField::setMaster()
{
    if (master_ != -1)
    {
        return;
    }

    const volScalarField& omega =
        static_cast<const volScalarField&>(this->internalField());

    const volScalarField::Boundary& bf = omega.boundaryField();

    label master = -1;
    forAll(bf, patchi)
    {
        if (isA<omegaRoughWallFunctionFvPatchScalarField>(bf[patchi]))
        {
            omegaRoughWallFunctionFvPatchScalarField& opf = omegaPatch(patchi);

            if (master == -1)
            {
                master = patchi;
            }

            opf.master() = master;
        }
    }
}


void omegaRoughWallFunctionFvPatchScalarField::createAveragingWeights()
{
    const volScalarField& omega =
        static_cast<const volScalarField&>(this->internalField());

    const volScalarField::Boundary& bf = omega.boundaryField();

    const fvMesh& mesh = omega.mesh();

    if (initialised_ && !mesh.changing())
    {
        return;
    }

    volScalarField weights
    (
        IOobject
        (
            "weights",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false // do not register
        ),
        mesh,
        dimensionedScalar(dimless, 0)
    );

    DynamicList<label> omegaPatches(bf.size());
    forAll(bf, patchi)
    {
        if (isA<omegaRoughWallFunctionFvPatchScalarField>(bf[patchi]))
        {
            omegaPatches.append(patchi);

            const labelUList& faceCells = bf[patchi].patch().faceCells();
            forAll(faceCells, i)
            {
                label celli = faceCells[i];
                weights[celli]++;
            }
        }
    }

    cornerWeights_.setSize(bf.size());
    forAll(omegaPatches, i)
    {
        label patchi = omegaPatches[i];
        const fvPatchScalarField& wf = weights.boundaryField()[patchi];
        cornerWeights_[patchi] = 1.0/wf.patchInternalField();
    }

    G_.setSize(internalField().size(), 0.0);
    omega_.setSize(internalField().size(), 0.0);

    initialised_ = true;
}


omegaRoughWallFunctionFvPatchScalarField&
omegaRoughWallFunctionFvPatchScalarField::omegaPatch(const label patchi)
{
    const volScalarField& omega =
        static_cast<const volScalarField&>(this->internalField());

    const volScalarField::Boundary& bf = omega.boundaryField();

    const omegaRoughWallFunctionFvPatchScalarField& opf =
        refCast<const omegaRoughWallFunctionFvPatchScalarField>(bf[patchi]);

    return const_cast<omegaRoughWallFunctionFvPatchScalarField&>(opf);
}


void omegaRoughWallFunctionFvPatchScalarField::calculateTurbulenceFields
(
    const momentumTransportModel& turbModel,
    scalarField& G0,
    scalarField& omega0
)
{
    // accumulate all of the G and omega contributions
    forAll(cornerWeights_, patchi)
    {
        if (!cornerWeights_[patchi].empty())
        {
            omegaRoughWallFunctionFvPatchScalarField& opf = omegaPatch(patchi);

            const List<scalar>& w = cornerWeights_[patchi];

            opf.calculate(turbModel, w, opf.patch(), G0, omega0);
        }
    }

    // apply zero-gradient condition for omega
    forAll(cornerWeights_, patchi)
    {
        if (!cornerWeights_[patchi].empty())
        {
            omegaRoughWallFunctionFvPatchScalarField& opf = omegaPatch(patchi);

            opf == scalarField(omega0, opf.patch().faceCells());
        }
    }
}


void omegaRoughWallFunctionFvPatchScalarField::calculate
(
    const momentumTransportModel& turbModel,
    const List<scalar>& cornerWeights,
    const fvPatch& patch,
    scalarField& G0,
    scalarField& omega0
)
{
    const label patchi = patch.index();

    const nutWallFunctionFvPatchScalarField& nutw =
        nutWallFunctionFvPatchScalarField::nutw(turbModel, patchi);

    const scalarField& y = turbModel.y()[patchi];

    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();

    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];

    const scalarField magGradUw(mag(Uw.snGrad()));

    const FieldType& G =
        db().lookupObject<FieldType>(turbModel.GName());

    const scalar Cmu25 = pow025(nutw.Cmu());
    const scalar Cmu5 = sqrt(nutw.Cmu());
    
    label nRoughCells = 0;
    scalar KsPlusAve = 0; //using for debug
    const scalar beta5 = pow(0.075, 0.5);
    const scalar betaOmega_ = 0.0708; //TODO define in .H
    

    // Set omega and G
    forAll(nutw, facei)
    {
        const label celli = patch.faceCells()[facei];
        const scalar w = cornerWeights[facei];

        const scalar Rey = y[facei]*sqrt(k[celli])/nuw[facei];
        const scalar yPlus = Cmu25*Rey;
        const scalar uPlus = (1/nutw.kappa())*log(nutw.E()*yPlus);

        if (blended_)
        {
            const scalar lamFrac = exp(-Rey/11);
            const scalar turbFrac = 1 - lamFrac;

            const scalar uStar = sqrt
            (
                lamFrac*nuw[facei]*magGradUw[facei] + turbFrac*Cmu5*k[celli]
            );

            const scalar omegaVis = 6*nuw[facei]/(beta1_*sqr(y[facei]));
            const scalar omegaLog = uStar/(Cmu5*nutw.kappa()*y[facei]);

            omega0[celli] += w*(lamFrac*omegaVis + turbFrac*omegaLog);

            G0[celli] +=
                w
               *(
                   lamFrac*G[celli]

                 + turbFrac
                  *sqr(uStar*magGradUw[facei]*y[facei]/uPlus)
                  /(nuw[facei]*nutw.kappa()*yPlus)
               );
        }
        else
        {
            //if (yPlus < nutw.yPlusLam())
            //{
            	//Adding roughness
            	const scalar ut = sqrt( (nuw[facei]+nutw[facei]) * magGradUw[facei] );
            	const scalar KsPlus = Ks_[facei]*ut / nuw[facei];
            	KsPlusAve+=KsPlus;
            	if (KsPlus > 5) // if KsPlus less then 5 then surface is considered smooth but this check is not needed            	
            	{
            		nRoughCells++;
		}
			
		    	const scalar ds = 0.03 * Ks_[facei] * min( 1., pow(KsPlus/30., 2./3.))
							     * min( 1., pow(KsPlus/45., 1./4.))
							     * min( 1., pow(KsPlus/60., 1./4.));
			
			const scalar omegaVis = min(
					    ut  / ( beta5 * nutw.kappa() * ds),
					    6.*nuw[facei]/betaOmega_/y[facei]/y[facei]
					    );							
            /*	}
            	else
            	{
		        const scalar omegaVis = 6*nuw[facei]/(beta1_*sqr(y[facei]));

                }
             */   
                omega0[celli] += w*omegaVis;
		 G0[celli] += w*G[celli];
            //}
	    /*
            else
            {
                const scalar uStar = sqrt(Cmu5*k[celli]);
                const scalar omegaLog = uStar/(Cmu5*nutw.kappa()*y[facei]);

                omega0[celli] += w*omegaLog;

                G0[celli] +=
                    w*
                    sqr(uStar*magGradUw[facei]*y[facei]/uPlus)
                   /(nuw[facei]*nutw.kappa()*yPlus);
            }
	    */
        }
    }
    
    if (true)
    {
        Pout<< "Patch: " << patch.name()
            << ": number of rough faces = " << nRoughCells
            << " out of " << patch.size()
            << " with average ks+ = " << KsPlusAve/patch.size()
            << endl;
    }
    
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

omegaRoughWallFunctionFvPatchScalarField::omegaRoughWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF),
    Ks_(p.size(), 0.0),
    beta1_(0.075),
    blended_(false),
    G_(),
    omega_(),
    initialised_(false),
    master_(-1),
    cornerWeights_()
{}


omegaRoughWallFunctionFvPatchScalarField::omegaRoughWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict),
    Ks_("Ks", dict, p.size()),
    beta1_(dict.lookupOrDefault<scalar>("beta1", 0.075)),
    blended_(dict.lookupOrDefault<Switch>("blended", false)),
    G_(),
    omega_(),
    initialised_(false),
    master_(-1),
    cornerWeights_()
{
    // apply zero-gradient condition on start-up
    this->operator==(patchInternalField());
}


omegaRoughWallFunctionFvPatchScalarField::omegaRoughWallFunctionFvPatchScalarField
(
    const omegaRoughWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper),
    Ks_(mapper(ptf.Ks_)),
    beta1_(ptf.beta1_),
    blended_(ptf.blended_),
    G_(),
    omega_(),
    initialised_(false),
    master_(-1),
    cornerWeights_()
{}


omegaRoughWallFunctionFvPatchScalarField::omegaRoughWallFunctionFvPatchScalarField
(
    const omegaRoughWallFunctionFvPatchScalarField& owfpsf
)
:
    fixedValueFvPatchField<scalar>(owfpsf),
    Ks_(owfpsf.Ks_),
    beta1_(owfpsf.beta1_),
    blended_(owfpsf.blended_),
    G_(),
    omega_(),
    initialised_(false),
    master_(-1),
    cornerWeights_()
{}


omegaRoughWallFunctionFvPatchScalarField::omegaRoughWallFunctionFvPatchScalarField
(
    const omegaRoughWallFunctionFvPatchScalarField& owfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(owfpsf, iF),
    Ks_(owfpsf.Ks_),
    beta1_(owfpsf.beta1_),
    blended_(owfpsf.blended_),
    G_(),
    omega_(),
    initialised_(false),
    master_(-1),
    cornerWeights_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalarField& omegaRoughWallFunctionFvPatchScalarField::G(bool init)
{
    if (patch().index() == master_)
    {
        if (init)
        {
            G_ = 0.0;
        }

        return G_;
    }

    return omegaPatch(master_).G();
}


scalarField& omegaRoughWallFunctionFvPatchScalarField::omega(bool init)
{
    if (patch().index() == master_)
    {
        if (init)
        {
            omega_ = 0.0;
        }

        return omega_;
    }

    return omegaPatch(master_).omega(init);
}


void omegaRoughWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const momentumTransportModel& turbModel =
        db().lookupObject<momentumTransportModel>
        (
            IOobject::groupName
            (
                momentumTransportModel::typeName,
                internalField().group()
            )
        );

    setMaster();

    if (patch().index() == master_)
    {
        createAveragingWeights();
        calculateTurbulenceFields(turbModel, G(true), omega(true));
    }

    const scalarField& G0 = this->G();
    const scalarField& omega0 = this->omega();

    FieldType& G =
        const_cast<FieldType&>
        (
            db().lookupObject<FieldType>(turbModel.GName())
        );

    FieldType& omega = const_cast<FieldType&>(internalField());

    forAll(*this, facei)
    {
        label celli = patch().faceCells()[facei];

        G[celli] = G0[celli];
        omega[celli] = omega0[celli];
    }

    fvPatchField<scalar>::updateCoeffs();
}


void omegaRoughWallFunctionFvPatchScalarField::updateWeightedCoeffs
(
    const scalarField& weights
)
{
    if (updated())
    {
        return;
    }

    const momentumTransportModel& turbModel =
        db().lookupObject<momentumTransportModel>
        (
            IOobject::groupName
            (
                momentumTransportModel::typeName,
                internalField().group()
            )
        );

    setMaster();

    if (patch().index() == master_)
    {
        createAveragingWeights();
        calculateTurbulenceFields(turbModel, G(true), omega(true));
    }

    const scalarField& G0 = this->G();
    const scalarField& omega0 = this->omega();

    FieldType& G =
        const_cast<FieldType&>
        (
            db().lookupObject<FieldType>(turbModel.GName())
        );

    FieldType& omega = const_cast<FieldType&>(internalField());

    scalarField& omegaf = *this;

    // only set the values if the weights are > tolerance
    forAll(weights, facei)
    {
        scalar w = weights[facei];

        if (w > tolerance_)
        {
            label celli = patch().faceCells()[facei];

            G[celli] = (1.0 - w)*G[celli] + w*G0[celli];
            omega[celli] = (1.0 - w)*omega[celli] + w*omega0[celli];
            omegaf[facei] = omega[celli];
        }
    }

    fvPatchField<scalar>::updateCoeffs();
}


void omegaRoughWallFunctionFvPatchScalarField::manipulateMatrix
(
    fvMatrix<scalar>& matrix
)
{
    if (manipulatedMatrix())
    {
        return;
    }

    matrix.setValues(patch().faceCells(), patchInternalField());

    fvPatchField<scalar>::manipulateMatrix(matrix);
}


void omegaRoughWallFunctionFvPatchScalarField::manipulateMatrix
(
    fvMatrix<scalar>& matrix,
    const Field<scalar>& weights
)
{
    if (manipulatedMatrix())
    {
        return;
    }

    DynamicList<label> constraintCells(weights.size());
    DynamicList<scalar> constraintomega(weights.size());
    const labelUList& faceCells = patch().faceCells();

    const DimensionedField<scalar, volMesh>& omega
        = internalField();

    label nConstrainedCells = 0;


    forAll(weights, facei)
    {
        // only set the values if the weights are > tolerance
        if (weights[facei] > tolerance_)
        {
            nConstrainedCells++;

            label celli = faceCells[facei];

            constraintCells.append(celli);
            constraintomega.append(omega[celli]);
        }
    }

    if (debug)
    {
        Pout<< "Patch: " << patch().name()
            << ": number of constrained cells = " << nConstrainedCells
            << " out of " << patch().size()
            << endl;
    }

    matrix.setValues
    (
        constraintCells,
        scalarField(constraintomega)
    );

    fvPatchField<scalar>::manipulateMatrix(matrix);
}


void omegaRoughWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    writeEntry(os, "Ks", Ks_);
    writeEntry(os, "beta1", beta1_);
    writeEntry(os, "blended", blended_);
    fixedValueFvPatchField<scalar>::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    omegaRoughWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
