/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2015 OpenFOAM Foundation
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

#include "twoPhaseCavitationMixture.H"

#include "gradientEnergyFvPatchScalarField.H"
#include "mixedEnergyFvPatchScalarField.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(twoPhaseCavitationMixture, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twoPhaseCavitationMixture::twoPhaseCavitationMixture
(
    const fvMesh& mesh
)
:
    psiThermo(mesh, word::null),
    twoPhaseMixture(mesh, *this),
    thermo1_(NULL),
    thermo2_(NULL)
{
    {
        volScalarField T1(IOobject::groupName("T", phase1Name()), T_);
        T1.write();
    }

    {
        volScalarField T2(IOobject::groupName("T", phase2Name()), T_);
        T2.write();
    }

    thermo1_ = rhoThermo::New(mesh, phase1Name());
    thermo2_ = rhoThermo::New(mesh, phase2Name());

    thermo1_->validate(phase1Name(), "e");
    thermo2_->validate(phase2Name(), "e");

    correct();


    // TODO Adapt this code for the two phase model
    //Initializing the cavitation model
//	forAllIter(PtrDictionary<phaseModel>, phases_, phasei)
//	{
//		if (phasei().name() == "water")
//		{
//			//Info << "phasei name: " << phasei().name() << "\n";
//			const phaseModel& alpha1 = phasei();
//			tmp<volScalarField> TMPrho1(phasei().thermo().rho());
//			const volScalarField& rho1(TMPrho1());
//
//			forAllIter(PtrDictionary<phaseModel>, phases_, phasej)
//			{
//				if (phasej().name() == "vapor")
//				{
//					//Info << "phasej name: " << phasej().name() << "\n";
//					const phaseModel& alpha2 = phasej();
//					tmp<volScalarField> TMPrho2(phasej().thermo().rho());
//					const volScalarField& rho2(TMPrho2());
//
//					cavitationModel_ = MultiphaseCavitation::New(U, 	phi,
//																		rho1,
//																		alpha1,
//																		alpha2);
//				}
//			}
//		}
//	}
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::twoPhaseCavitationMixture::~twoPhaseCavitationMixture()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::twoPhaseCavitationMixture::correct()
{
    thermo1_->he() = thermo1_->he(p_, T_);
    thermo1_->correct();

    thermo2_->he() = thermo2_->he(p_, T_);
    thermo2_->correct();

    psi_ = alpha1()*thermo1_->psi() + alpha2()*thermo2_->psi();
    mu_ = alpha1()*thermo1_->mu() + alpha2()*thermo2_->mu();
    alpha_ = alpha1()*thermo1_->alpha() + alpha2()*thermo2_->alpha();
}


bool Foam::twoPhaseCavitationMixture::incompressible() const
{
    return thermo1_->incompressible() && thermo2_->incompressible();
}


bool Foam::twoPhaseCavitationMixture::isochoric() const
{
    return thermo1_->isochoric() && thermo2_->isochoric();
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseCavitationMixture::he
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    return alpha1()*thermo1_->he(p, T) + alpha2()*thermo2_->he(p, T);
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseCavitationMixture::he
(
    const scalarField& p,
    const scalarField& T,
    const labelList& cells
) const
{
    return
        scalarField(alpha1(), cells)*thermo1_->he(p, T, cells)
      + scalarField(alpha2(), cells)*thermo2_->he(p, T, cells);
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseCavitationMixture::he
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    return
        alpha1().boundaryField()[patchi]*thermo1_->he(p, T, patchi)
      + alpha2().boundaryField()[patchi]*thermo2_->he(p, T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseCavitationMixture::hc() const
{
    return alpha1()*thermo1_->hc() + alpha2()*thermo2_->hc();
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseCavitationMixture::THE
(
    const scalarField& h,
    const scalarField& p,
    const scalarField& T0,
    const labelList& cells
) const
{
    NotImplemented;
    return T0;
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseCavitationMixture::THE
(
    const scalarField& h,
    const scalarField& p,
    const scalarField& T0,
    const label patchi
) const
{
    NotImplemented;
    return T0;
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseCavitationMixture::Cp() const
{
    return alpha1()*thermo1_->Cp() + alpha2()*thermo2_->Cp();
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseCavitationMixture::Cp
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    return
        alpha1().boundaryField()[patchi]*thermo1_->Cp(p, T, patchi)
      + alpha2().boundaryField()[patchi]*thermo2_->Cp(p, T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseCavitationMixture::Cv() const
{
    return alpha1()*thermo1_->Cv() + alpha2()*thermo2_->Cv();
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseCavitationMixture::Cv
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    return
        alpha1().boundaryField()[patchi]*thermo1_->Cv(p, T, patchi)
      + alpha2().boundaryField()[patchi]*thermo2_->Cv(p, T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseCavitationMixture::gamma() const
{
    return alpha1()*thermo1_->gamma() + alpha2()*thermo2_->gamma();
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseCavitationMixture::gamma
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    return
        alpha1().boundaryField()[patchi]*thermo1_->gamma(p, T, patchi)
      + alpha2().boundaryField()[patchi]*thermo2_->gamma(p, T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseCavitationMixture::Cpv() const
{
    return alpha1()*thermo1_->Cpv() + alpha2()*thermo2_->Cpv();
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseCavitationMixture::Cpv
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    return
        alpha1().boundaryField()[patchi]*thermo1_->Cpv(p, T, patchi)
      + alpha2().boundaryField()[patchi]*thermo2_->Cpv(p, T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseCavitationMixture::CpByCpv() const
{
    return
        alpha1()*thermo1_->CpByCpv()
      + alpha2()*thermo2_->CpByCpv();
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseCavitationMixture::CpByCpv
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    return
        alpha1().boundaryField()[patchi]*thermo1_->CpByCpv(p, T, patchi)
      + alpha2().boundaryField()[patchi]*thermo2_->CpByCpv(p, T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseCavitationMixture::nu() const
{
    return mu()/(alpha1()*thermo1_->rho() + alpha2()*thermo2_->rho());
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseCavitationMixture::nu
(
    const label patchi
) const
{
    return
        mu(patchi)
       /(
            alpha1().boundaryField()[patchi]*thermo1_->rho(patchi)
          + alpha2().boundaryField()[patchi]*thermo2_->rho(patchi)
        );
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseCavitationMixture::kappa() const
{
    return alpha1()*thermo1_->kappa() + alpha2()*thermo2_->kappa();
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseCavitationMixture::kappa
(
    const label patchi
) const
{
    return
        alpha1().boundaryField()[patchi]*thermo1_->kappa(patchi)
      + alpha2().boundaryField()[patchi]*thermo2_->kappa(patchi);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseCavitationMixture::kappaEff
(
    const volScalarField& alphat
) const
{
    return
        alpha1()*thermo1_->kappaEff(alphat)
      + alpha2()*thermo2_->kappaEff(alphat);
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseCavitationMixture::kappaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return
        alpha1().boundaryField()[patchi]*thermo1_->kappaEff(alphat, patchi)
      + alpha2().boundaryField()[patchi]*thermo2_->kappaEff(alphat, patchi)
    ;
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseCavitationMixture::alphaEff
(
    const volScalarField& alphat
) const
{
    return
        alpha1()*thermo1_->alphaEff(alphat)
      + alpha2()*thermo2_->alphaEff(alphat);
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseCavitationMixture::alphaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return
        alpha1().boundaryField()[patchi]*thermo1_->alphaEff(alphat, patchi)
      + alpha2().boundaryField()[patchi]*thermo2_->alphaEff(alphat, patchi)
    ;
}


// ************************************************************************* //
