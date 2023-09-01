/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "interfaceProperties.H"
#include "alphaContactAngleTwoPhaseFvPatchScalarField.H"
#include "mathematicalConstants.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "fvcReconstruct.H"
#include "unitConversion.H"

#include "pointMesh.H"
#include "pointFields.H"
//#include "volPointInterpolation.H"
//#include "interpolatePointToCell.H"
#include "primitivePatchInterpolation.H"
#include "wallFvPatch.H"

#include "fvCFD.H"
#include "fixedValuePointPatchFields.H"
#include "fvPatchFields.H"
#include "fvsPatchFields.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Correction for the boundary condition on the unit normal nHat on
// walls to produce the correct contact angle.

// The dynamic contact angle is calculated from the component of the
// velocity on the direction of the interface, parallel to the wall.

void Foam::interfaceProperties::correctContactAngle
(
    surfaceVectorField::Boundary& nHatb,
    const surfaceVectorField::Boundary& gradAlphaf,
	volScalarField& thetaRad_
) const
{
	Info << "Contact angle correction." << endl;
	
    const fvMesh& mesh = alpha1_.mesh();
    const volScalarField::Boundary& abf = alpha1_.boundaryField();
	
    const fvBoundaryMesh& boundary = mesh.boundary();
    forAll(boundary, patchi)
    {
        if (isA<alphaContactAngleTwoPhaseFvPatchScalarField>(abf[patchi]))
        {
            alphaContactAngleTwoPhaseFvPatchScalarField& acap =
                const_cast<alphaContactAngleTwoPhaseFvPatchScalarField&>
                (
                    refCast<const alphaContactAngleTwoPhaseFvPatchScalarField>
                    (
                        abf[patchi]
                    )
                );
			
            fvsPatchVectorField& nHatp = nHatb[patchi];
			
			scalarField theta(degToRad() * acap.theta(U_.boundaryField()[patchi], nHatp));
			
			const fvPatch& curPatch = boundary[patchi];
			
			forAll(curPatch,facei)
			{	
				label faceCelli = curPatch.faceCells()[facei];	
				thetaRad_[faceCelli] =  theta[facei];
			}	
			
            const vectorField nf(boundary[patchi].nf());

            const scalarField a12(nHatp & nf);
            const scalarField b1(cos(theta));
			
            scalarField b2(nHatp.size());
            forAll(b2, facei)
            {
                b2[facei] = cos(acos(a12[facei]) - theta[facei]);
            }

            const scalarField det(1.0 - a12*a12);

            scalarField a((b1 - a12*b2)/(det+VSMALL));
            scalarField b((b2 - a12*b1)/(det+VSMALL));

            nHatp = a*nf + b*nHatp;
            nHatp /= (mag(nHatp) + deltaN_.value());
			
            acap.gradient() = (nf & nHatp)*mag(gradAlphaf[patchi]);
            acap.evaluate();
				
        }
    }
	thetaRad_.correctBoundaryConditions();
}

void Foam::interfaceProperties::smoothen //https://github.com/floquation/OF-kva_interfaceProperties/blob/master/curvatureModel/vofsmooth/vofsmooth.C
(
    volScalarField& smooth_func
) const
{
    const fvMesh& mesh = smooth_func.mesh();
    const surfaceVectorField& Sf = mesh.Sf();

    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();
	int numSmoothingIterations_= 1;

    for(int iter = 0; iter < numSmoothingIterations_; iter++)
    {
    	scalarField smooth_cal(mesh.nCells(),scalar(0));

    	scalarField sum_area(mesh.nCells(),scalar(0));

		surfaceScalarField smoothF = fvc::interpolate(smooth_func);

		for(int facei = 0; facei < nei.size(); facei++) //KVA note: should be own???
		{
			smooth_cal[own[facei]] += smoothF[facei]*mag(Sf[facei]);
			sum_area[own[facei]] += mag(Sf[facei]);
		}

		forAll(nei,facei)
		{
			smooth_cal[nei[facei]] += smoothF[facei]*mag(Sf[facei]);
			sum_area[nei[facei]] += mag(Sf[facei]);
		}

		forAll(mesh.boundary(), patchi)
		{
			
			const labelList& pFaceCells = mesh.boundary()[patchi].faceCells(); //const unallocLabelList& pFaceCells = mesh.boundary()[patchi].faceCells();

			const fvsPatchScalarField& pssf = smoothF.boundaryField()[patchi];

			forAll(mesh.boundary()[patchi], facei)
			{
			   smooth_cal[pFaceCells[facei]] += pssf[facei]*mag(Sf[facei]);
			   sum_area[pFaceCells[facei]] += mag(Sf[facei]);
			}
		}

		forAll(mesh.cells(),celli)
		{
			smooth_func[celli] = smooth_cal[celli]/sum_area[celli];
		}

		smooth_func.correctBoundaryConditions();
    }
}

void Foam::interfaceProperties::calculateHeaviside()
{
	// Heaviside
	Info<<"   Heaviside"<<endl;
	
	forAll(psi_,celli)//forAll(psi_.mesh().cells(),celli)
    {
		if(psi_[celli] < -epsilon_.value())
			Heaviside_[celli] = double(0.0);
		else if(psi_[celli] > epsilon_.value())
			Heaviside_[celli] = double(1.0);
		else
			Heaviside_[celli] = double(1.0)/double(2.0)*(double(1.0)+psi_[celli]/epsilon_.value()+sin(M_PI*psi_[celli]/epsilon_.value())/M_PI);
    }
	Heaviside_.correctBoundaryConditions();
	
	const surfaceScalarField psif(fvc::interpolate(psi_));
	forAll(psif,facei)//forAll(psif.mesh().faces(),facei)
    {
		if(psif[facei] < -epsilon_.value())
			Heavisidef_[facei] = double(0.0);
		else if(psif[facei] > epsilon_.value())
			Heavisidef_[facei] = double(1.0);
		else
			Heavisidef_[facei] = double(1.0)/double(2.0)*(double(1.0)+psif[facei]/epsilon_.value()+sin(M_PI*psif[facei]/epsilon_.value())/M_PI);
    }
	
}


void Foam::interfaceProperties::calculateDirac()
{
	Info<<"   Dirac "<<endl;
	
	forAll(psi_,celli)//forAll(psi_.mesh().cells(),celli)
    {  
		if(mag(psi_[celli]) > epsilon_.value())
			Dirac_[celli] = double(0.0);
		else
			Dirac_[celli] = double(1.0)/(double(2.0)*epsilon_.value())*(double(1.0)+cos(M_PI*psi_[celli]/epsilon_.value()));
    }
	Dirac_.correctBoundaryConditions();
	
	const surfaceScalarField psif(fvc::interpolate(psi_));
	forAll(psif,facei)//forAll(psif.mesh().faces(),facei)
    {  
		if(mag(psif[facei]) > epsilon_.value())
			Diracf_[facei] = double(0.0);
		else
			Diracf_[facei] = double(1.0)/(double(2.0)*epsilon_.value())*(double(1.0)+cos(M_PI*psif[facei]/epsilon_.value()));
    }
	
}

void Foam::interfaceProperties::calculateK()
{
	Info << "Interface curvature" << endl;
		
	
	// Mesh 
    const fvMesh& mesh = alpha1_.mesh();
    const surfaceVectorField& Sf = mesh.Sf();
	const surfaceScalarField& magSf = mesh.magSf();
	
	const polyBoundaryMesh& boundary = mesh.boundaryMesh();

	//n from alpha1
	Info << "	n from alpha1" << endl;
	const volVectorField gradAlpha(fvc::grad(alpha1_, "nHat"));
	//volVectorField nHatAlpha_(gradAlpha/(mag(gradAlpha) + deltaN_));
	
	surfaceVectorField gradAlphaf(fvc::interpolate(gradAlpha));
	nHatAlphafv_ = gradAlphaf/(mag(gradAlphaf) + deltaN_);
	
	
	//n from Heaviside
	Info << "	n from Heaviside" << endl;
	const volVectorField gradHeaviside(fvc::grad(Heaviside_, "nHat"));
	//volVectorField nHatHeaviside_(gradHeaviside/(mag(gradHeaviside) + deltaN_));
	
	surfaceVectorField gradHeavisidef(fvc::interpolate(gradHeaviside));
	nHatHeavisidefv_ = gradHeavisidef/(mag(gradHeavisidef) + deltaN_);
	
	
	//n from psi
	Info << "	n from psi" << endl;
	const volVectorField gradPsi(fvc::grad(psi_, "nHat"));
	volVectorField nHatPsi_(gradPsi/(mag(gradPsi) + SMALL));

	const volVectorField gradPsiW(gradPsi*((Dirac_*epsilon_)+SMALL)/(mag(gradPsi)+SMALL));
	const surfaceVectorField gradPsiWf(fvc::interpolate(gradPsiW));
	surfaceVectorField nHatPsifv_(gradPsiWf/(mag(gradPsiWf) + SMALL));
	
	/////////SMOOTH !!!!
	/*
	const volVectorField gradPsiW(gradPsi*((Dirac_*epsilon_)+SMALL)/(mag(gradPsi)+SMALL));
	volVectorField gradPsiWSmooth(gradPsiW);
	
	for(int iter = 0; iter < 2; iter++)
    {
		Info << "	n smoothing " << iter+1 <<endl;
			
		volScalarField gradPsiWSmoothX = gradPsiWSmooth.component(0);
		volScalarField gradPsiWSmoothY = gradPsiWSmooth.component(1);
		volScalarField gradPsiWSmoothZ = gradPsiWSmooth.component(2);
		smoothen(gradPsiWSmoothX);
		smoothen(gradPsiWSmoothY);
		smoothen(gradPsiWSmoothZ);
		
		//gradPsiWSmooth.component(0) = gradPsiWSmoothX;
		//gradPsiWSmooth.component(1) = gradPsiWSmoothY;
		//gradPsiWSmooth.component(2) = gradPsiWSmoothZ;
		forAll(gradPsiWSmooth, cellI)
		{
			gradPsiWSmooth[cellI][0] = gradPsiWSmoothX[cellI];
			gradPsiWSmooth[cellI][1] = gradPsiWSmoothY[cellI];
			gradPsiWSmooth[cellI][2] = gradPsiWSmoothZ[cellI];
		}
		gradPsiWSmooth.correctBoundaryConditions();
		
		volScalarField gradPsiWSmoothMag = mag(gradPsiWSmooth);
		volScalarField gradPsiWMag = mag(gradPsiW);
		
		Info << gMax(gradPsiWSmoothMag) - gMax(gradPsiWMag) <<endl;

		
	}
	nHatPsi_ = gradPsiWSmooth/(mag(gradPsiWSmooth)+SMALL);
	nHatPsi_.correctBoundaryConditions();
	const surfaceVectorField gradPsiWf(fvc::interpolate(gradPsiWSmooth));
	surfaceVectorField nHatPsifv_(gradPsiWf/(mag(gradPsiWf) + SMALL));
	*/
	
	Info << "	Interface normal" << endl;
	nHat_ == nHatPsi_;
	nHatfv_ = nHatPsifv_;
	
	// coefficient for the compression strenght
	//compStrCoeff_ = (cos(double(2)*acos(mag(nHatAlphafv_ & (Sf/magSf))))+double(1.0))/double(2);
	/*
	compStrCoeff_ = (cos(double(2)*acos(mag(nHatHeavisidefv_ & (Sf/magSf))))+double(1.0))/double(2);
	
	forAll(compStrCoeff_,facei)
	{
		if(compStrCoeff_[facei]>double(1))
		{
			compStrCoeff_[facei]=double(1);
		}
	}
	*/
	const surfaceScalarField psif(fvc::interpolate(psi_));
	compStrCoeff_ = (cos(double(2.0)*acos(mag(nHatPsifv_ & (Sf/magSf))))+double(1.0))/double(2.0);
	forAll(compStrCoeff_,facei)
	{
		if(compStrCoeff_[facei]>double(1.0))
		{
			compStrCoeff_[facei]=double(1.0);
		}
		
		if(psif[facei]/epsilon_.value()>double(1.0))
		{
			compStrCoeff_[facei]=double(1.0);
		}
	}
		
	// Correct CA
	Info<<"   Contact angle"<<endl;
	thetaRad_ = thetaRad_*double(0) + constant::mathematical::pi/2;
	//correctContactAngle(nHatfv_.boundaryFieldRef(), gradAlphaf.boundaryField(),thetaRad_);
	correctContactAngle(nHatAlphafv_.boundaryFieldRef(), gradAlphaf.boundaryField(),thetaRad_);
	
	//Curvature
	//volScalarField Heaviside_smooth = Heaviside_;
	//smoothen(Heaviside_smooth);
	//const volVectorField gradHeavisideSmooth(fvc::grad(Heaviside_smooth, "nHat"));
	//surfaceVectorField gradHeavisideSmoothf(fvc::interpolate(gradHeavisideSmooth));
	//surfaceVectorField nHatSmoothfv_ = gradHeavisideSmoothf/(mag(gradHeavisideSmoothf) + deltaN_);
	//surfaceScalarField nHatf_(nHatSmoothfv_ & Sf); 
	
	// Face unit interface normal flux
    
	surfaceScalarField nHatf_(nHatfv_ & Sf);
		
    // Simple expression for curvature
    K_ = -fvc::div(nHatf_);
	
	/*
    // Complex expression for curvature.
    // Correction is formally zero but numerically non-zero.
    
    forAll(nHat_.boundaryField(), patchi)
    {
        nHat_.boundaryFieldRef()[patchi] = nHatfv_.boundaryField()[patchi];
    }

    K_ = -fvc::div(nHatf_) + (nHat_ & fvc::grad(nHatfv_) & nHat_);
	*/
}	


void Foam::interfaceProperties::calculateSurfaceTension()
{
	sigmaIP_ = sigmaPtr_->sigma();
}	
	


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceProperties::interfaceProperties
(
    const volScalarField& alpha1,
    const volVectorField& U,
	const volScalarField& psi,
	const dimensionedScalar& epsilon,
    const IOdictionary& dict
)
:
    transportPropertiesDict_(dict),
	
    cAlpha_
    (
        alpha1.mesh().solverDict(alpha1.name()).get<scalar>("cAlpha")
    ),
	
    sigmaPtr_(surfaceTensionModel::New(dict, alpha1.mesh())),
	
	epsilon_(epsilon),
	
    deltaN_
    (
        "deltaN",
        1e-8/epsilon_//1e-8/cbrt(average(alpha1.mesh().V())) // Mesh width min
    ),
	
    alpha1_(alpha1),
	
    U_(U),
	
	psi_(psi),	
	
	// Interface normal
	
	nHat_
    (
        IOobject
        (
            "nHat",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        fvc::grad(alpha1_)/(mag(fvc::grad(alpha1_))+deltaN_)
    ),
	
	nHatfv_
	(
        IOobject
        (
            "nHatfv",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        fvc::interpolate(nHat_)
    ),
	
	nHatAlphafv_
	(
        IOobject
        (
            "nHatAlphafv",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        fvc::interpolate(nHat_)
    ),
	
	nHatHeavisidefv_
	(
        IOobject
        (
            "nHatHeavisidefv_",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        fvc::interpolate(nHat_)
    ),
		
	compStrCoeff_
	(
        IOobject
        (
            "compStrCoeff",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
		fvc::interpolate(alpha1_)
    ),
	
	Heaviside_
	(
        IOobject
        (
            "Heaviside",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
		alpha1_
    ),
	
	Heavisidef_
	(
        IOobject
        (
            "Heavisidef",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
		fvc::interpolate(Heaviside_)
    ),
	
	Dirac_
	(
        IOobject
        (
            "Dirac",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
		mag(fvc::grad(Heaviside_))
    ),
	
	Diracf_
	(
        IOobject
        (
            "Diracf",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
		fvc::interpolate(Dirac_)
    ),
	
	K_
    (
        IOobject
        (
            "interfaceProperties:K",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_.mesh(),
        dimensionedScalar(dimless/dimLength, Zero)
    ),
	
	thetaRad_
	(
        IOobject
        (
            "thetaRad",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_.mesh(),
        dimensionedScalar(dimless, Zero)
    ),
	
	sigmaIP_
    (
        IOobject
        (
            "sigmaIP",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_.mesh(),
        dimensionedScalar(dimMass/(dimTime*dimTime), Zero)
    )
{
	calculateHeaviside();
	calculateDirac();
	calculateK();
	calculateSurfaceTension();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::interfaceProperties::sigmaK() const
{
	Info << "sigmaK calculation" << endl; 
    return sigmaIP_*K_;
}

Foam::tmp<Foam::volVectorField>
Foam::interfaceProperties::surfaceTensionForce() const
{
	Info << "Surface tension force calculation (volVectorField)" << endl;
	return sigmaK()*Dirac_*nHat_;
}

Foam::tmp<Foam::volVectorField>
Foam::interfaceProperties::surfaceTensionDensityScaledForce() const
{
	Info << "Surface tension force calculation (volVectorField)" << endl;
	return double(2)*Heaviside_ *sigmaK()*Dirac_*nHat_;
}

Foam::tmp<Foam::surfaceScalarField>
Foam::interfaceProperties::surfaceTensionForcef() const
{
	Info << "Surface tension force calculation (surfaceScalarField)" << endl;
	const fvMesh& mesh = alpha1_.mesh();
	return fvc::interpolate(sigmaK()*Dirac_)*(nHatfv_ & (mesh.Sf()/mesh.magSf()));
}

Foam::tmp<Foam::surfaceScalarField>
Foam::interfaceProperties::surfaceTensionDensityScaledForcef() const
{
	Info << "Surface tension force calculation (surfaceScalarField)" << endl;
	const fvMesh& mesh = alpha1_.mesh();
	return fvc::interpolate(double(2)*Heaviside_*sigmaK()*Dirac_)*(nHatfv_ & (mesh.Sf()/mesh.magSf()));
	//return double(2)*Heavisidef_*Diracf_*fvc::interpolate(sigmaK())*(nHatfv_ & (mesh.Sf()/mesh.magSf()));
}

Foam::tmp<Foam::volVectorField>
Foam::interfaceProperties::MarangoniForce() const
{
	volVectorField gradSigma(fvc::grad(sigmaIP_));
	volVectorField gradSigmaTan = gradSigma - nHat_*(nHat_ & gradSigma);
	
	return Dirac_ * gradSigmaTan;
}

Foam::tmp<Foam::volVectorField>
Foam::interfaceProperties::MarangoniDensityScaledForce() const
{
	volVectorField gradSigma(fvc::grad(sigmaIP_));
	volVectorField gradSigmaTan = gradSigma - nHat_*(nHat_ & gradSigma);
	
	return double(2)*Heaviside_ * Dirac_ * gradSigmaTan;
}

Foam::tmp<Foam::surfaceScalarField>
Foam::interfaceProperties::MarangoniForcef() const
{
	Info << "Marangoni force calculation (surfaceScalarField)" << endl;
	
	const fvMesh& mesh = alpha1_.mesh();
	volVectorField gradSigma(fvc::grad(sigmaIP_));
	surfaceVectorField gradSigmaf(fvc::interpolate(gradSigma));
	surfaceVectorField gradSigmaTanf = gradSigmaf - nHatfv_*(nHatfv_ & gradSigmaf);
	
    return fvc::interpolate(Dirac_)* (gradSigmaTanf & (mesh.Sf()/mesh.magSf()));
}


Foam::tmp<Foam::surfaceScalarField>
Foam::interfaceProperties::MarangoniDensityScaledForcef() const
{
	Info << "Marangoni force calculation (surfaceScalarField)" << endl;
	
	const fvMesh& mesh = alpha1_.mesh();
	volVectorField gradSigma(fvc::grad(sigmaIP_));
	surfaceVectorField gradSigmaf(fvc::interpolate(gradSigma));
	surfaceVectorField gradSigmaTanf = gradSigmaf - nHatfv_*(nHatfv_ & gradSigmaf);
	
    return fvc::interpolate(double(2)*Heaviside_*Dirac_)* (gradSigmaTanf & (mesh.Sf()/mesh.magSf()));
	//return double(2)*Heavisidef_*Diracf_* (gradSigmaTanf & (mesh.Sf()/mesh.magSf()));
}

Foam::tmp<Foam::volScalarField>
Foam::interfaceProperties::nearInterface() const
{
    return pos0(alpha1_ - 0.01)*pos0(0.99 - alpha1_);
}

Foam::tmp<Foam::volScalarField>
Foam::interfaceProperties::nearInterface1Epsilon() const
{
	return pos0((epsilon_ - mag(psi_))/epsilon_);
}


void Foam::interfaceProperties::correct()
{
	calculateHeaviside();
	calculateDirac();
	calculateK();
	calculateSurfaceTension();
}

bool Foam::interfaceProperties::read()
{
	
    alpha1_.mesh().solverDict(alpha1_.name()).readEntry("cAlpha", cAlpha_);
    sigmaPtr_->readDict(transportPropertiesDict_);
	
    return true;
}


// ************************************************************************* //
