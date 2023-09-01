/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.
    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.
    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.
    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
\*---------------------------------------------------------------------------*/

#include "saaDynamicKistlerAlphaContactAngleFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "fvPatchFields.H"
#include "volMesh.H"
#include "surfaceInterpolate.H"
#include "fvcReconstruct.H" 

//https://reader.elsevier.com/reader/sd/pii/S0001868614002280?token=E1DAF9CA630C6DAF7BC8BFF3F99FD0D89E35F40A732CE263D61F258D430257008F0A722329A349C0B7DAB0CCF1CD07DF&originRegion=us-east-1&originCreation=20220428010204
// https://riunet.upv.es/bitstream/handle/10251/99852/5020-15575-1-PB.pdf?sequence=1&isAllowed=y
// https://www.sciencedirect.com/science/article/pii/S036031991830394X

// * * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * //
namespace Foam
{
const scalar saaDynamicKistlerAlphaContactAngleFvPatchScalarField::convertToDeg =
    180.0/constant::mathematical::pi;

const scalar saaDynamicKistlerAlphaContactAngleFvPatchScalarField::convertToRad =
    constant::mathematical::pi/180.0;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::saaDynamicKistlerAlphaContactAngleFvPatchScalarField::
saaDynamicKistlerAlphaContactAngleFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(p, iF),
    thetaA_(0.0),
    thetaR_(0.0)
{}

Foam::saaDynamicKistlerAlphaContactAngleFvPatchScalarField::
saaDynamicKistlerAlphaContactAngleFvPatchScalarField
(
    const saaDynamicKistlerAlphaContactAngleFvPatchScalarField& gcpsf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(gcpsf, p, iF, mapper),
    thetaA_(gcpsf.thetaA_),
    thetaR_(gcpsf.thetaR_)
{}


Foam::saaDynamicKistlerAlphaContactAngleFvPatchScalarField::
saaDynamicKistlerAlphaContactAngleFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(p, iF),
    thetaA_(readScalar(dict.lookup("thetaA"))),
    thetaR_(readScalar(dict.lookup("thetaR")))
{
    evaluate();
}


Foam::saaDynamicKistlerAlphaContactAngleFvPatchScalarField::
saaDynamicKistlerAlphaContactAngleFvPatchScalarField
(
    const saaDynamicKistlerAlphaContactAngleFvPatchScalarField& gcpsf
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(gcpsf),
    thetaA_(gcpsf.thetaA_),
    thetaR_(gcpsf.thetaR_)
{}


Foam::saaDynamicKistlerAlphaContactAngleFvPatchScalarField::
saaDynamicKistlerAlphaContactAngleFvPatchScalarField
(
    const saaDynamicKistlerAlphaContactAngleFvPatchScalarField& gcpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(gcpsf, iF),
    thetaA_(gcpsf.thetaA_),
    thetaR_(gcpsf.thetaR_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> 
Foam::saaDynamicKistlerAlphaContactAngleFvPatchScalarField::theta
(
    const fvPatchVectorField& Up,
    const fvsPatchVectorField& nHat
) const
{
	
	Info << "Kistler model" << endl;
	
	
	const dictionary& transportProperties = db().lookupObject<IOdictionary>("transportProperties");
	scalar nu1(readScalar(transportProperties.subDict("water").lookup("nu")));
	scalar rho1(readScalar(transportProperties.subDict("water").lookup("rho")));
	scalar mup = nu1* rho1;
	
	scalar thetaMin_(transportProperties.subDict("water").lookupOrDefault<double>("thetaMin", 10));
	scalar thetaMax_(transportProperties.subDict("water").lookupOrDefault<double>("thetaMax", 170));
	
	//thetaMin_ = max(thetaMin_,SMALL);
	//thetaMax_ = max(thetaMax_,180-SMALL);
	scalar thetaMin(max(thetaMin_*convertToRad,SMALL));
	scalar thetaMax(min(thetaMax_*convertToRad,constant::mathematical::pi-SMALL));
	
	//scalar thetaMin = scalar(0.0)+convertToRad*SMALL; //convertToRad*double(45);
	//scalar thetaMax = constant::mathematical::pi-convertToRad*SMALL; //convertToRad*double(170);
	
	
	
	//Info << db().names() << endl;
	
	//const fvPatchField<scalar>& mup = patch().lookupPatchField<volScalarField, scalar>("mu");
	//const fvPatchField<scalar>& sigmap = patch().lookupPatchField<volScalarField, scalar>("sigmaIP");
	//const fvPatchField<vector>& UpSmoothp = patch().lookupPatchField<volVectorField, vector>("UKistler");
	//const fvsPatchField<vector>& nHatKistlerp = patch().lookupPatchField<surfaceVectorField, vector>("nHatKistler");
	
	//const polyPatch& pp = alphap.patch().patch();
	
	
	const fvPatchField<scalar>& sigmap = patch().lookupPatchField<volScalarField, scalar>("sigmaIP");

	//const fvPatchField<scalar>& wallCACorrp = patch().lookupPatchField<volScalarField, scalar>("wallCACorr");

    const vectorField nf(patch().nf());

    // Calculate the component of the velocity parallel to the wall
	vectorField Uwall(Up.patchInternalField() - Up);
    Uwall -= (nf & Uwall)*nf;

    // Find the direction of the interface parallel to the wall
	vectorField nWall(nHat - (nf & nHat)*nf);
	
    // Normalise nWall
	nWall /= (mag(nWall) + SMALL);
	
    // Calculate Uwall resolved normal to the interface parallel to
    // the wall
	scalarField uwall(nWall & Uwall);

	Info << "Ca number " << endl;
    //eb - Calculate local Capillary number
    scalarField Ca(-mup*uwall/ (sigmap + SMALL)); //scalarField Ca = mup*mag(uwall)/sigmap;

    //eb - Instantiate function object InverseHoffmanFunction for thetaA and thetaR
	
	scalar coefA = AngleCoefficient(convertToRad*thetaA_);
	scalar coefR = AngleCoefficient(convertToRad*thetaR_);
	scalar thetaS = acos((coefA*cos(convertToRad*thetaA_) + coefR*cos(convertToRad*thetaR_))/(coefA+coefR));
	
	/*
    Foam::saaDynamicKistlerAlphaContactAngleFvPatchScalarField::InverseHoffmanFunction
    InvHoffFuncThetaA
    (
        convertToRad*thetaA_
    );
    Foam::saaDynamicKistlerAlphaContactAngleFvPatchScalarField::InverseHoffmanFunction
    InvHoffFuncThetaR
    (
        convertToRad*thetaR_
    );
	*/
	Foam::saaDynamicKistlerAlphaContactAngleFvPatchScalarField::InverseHoffmanFunction
    InvHoffFuncThetaS
    (
        thetaS
    );
	Foam::saaDynamicKistlerAlphaContactAngleFvPatchScalarField::InverseHoffmanFunction
    InvHoffFuncThetaSMALL
    (
        convertToRad*SMALL
    );

    //eb - Calculate InverseHoffmanFunction for thetaA and thetaR using
    // RiddersRoot
	/*
    RiddersRoot RRInvHoffFuncThetaA(InvHoffFuncThetaA, 1.e-10);
    scalar InvHoffFuncThetaAroot = RRInvHoffFuncThetaA.root(0,65);
    RiddersRoot RRInvHoffFuncThetaR(InvHoffFuncThetaR, 1.e-10);
    scalar InvHoffFuncThetaRroot = RRInvHoffFuncThetaR.root(0,65);
	*/
	RiddersRoot RRInvHoffFuncThetaS(InvHoffFuncThetaS, 1.e-10);
    scalar InvHoffFuncThetaSroot = RRInvHoffFuncThetaS.root(0,65);
	RiddersRoot RRInvHoffFuncThetaSMALL(InvHoffFuncThetaSMALL, 1.e-10);
    scalar InvHoffFuncThetaSMALLroot = RRInvHoffFuncThetaSMALL.root(0,65);

    //eb - Calculate and return the value of contact angle on patch faces, a general approach: the product of Uwall and nWall is negative for advancing and positiv for receding motion.
    //     thetaDp is initialized to the equilibrium contact angle https://pubs.acs.org/doi/10.1021/la049410h
	
	Info << "Equilibrium CA: " << convertToDeg*thetaS << " [deg]" << endl;
    scalarField thetaDp(patch().size(), constant::mathematical::pi/2);
    forAll(uwall, pfacei)
    {
		scalar thetaK = constant::mathematical::pi/2;
		//if(mag(psip[pfacei])<(1-SMALL)*double(2)*epsilonp[pfacei])
		//if(wallCACorrp[pfacei]>SMALL)
		//{
			if(uwall[pfacei] < 0.0)
			{
				thetaK = max(HoffmanFunction(Ca[pfacei] + InvHoffFuncThetaSroot),convertToRad*thetaA_);
			}
			else if (uwall[pfacei] > 0.0)
			{
				thetaK = min(HoffmanFunction(max(Ca[pfacei] + InvHoffFuncThetaSroot,InvHoffFuncThetaSMALLroot)),convertToRad*thetaR_);
			}
		//}
		
		//thetaDp[pfacei] = min(max(thetaK,scalar(0.0)+convertToRad*SMALL),constant::mathematical::pi-convertToRad*SMALL);
		thetaDp[pfacei] = min(max(thetaK,thetaMin),thetaMax);
	}
	

    
    return convertToDeg*thetaDp;
}


Foam::scalar Foam::saaDynamicKistlerAlphaContactAngleFvPatchScalarField::HoffmanFunction
(
    const scalar& x
) const
{
    return acos(1 - 2*tanh(5.16*pow(x/(1+1.31*pow(x,0.99)),0.706)));
}

Foam::scalar Foam::saaDynamicKistlerAlphaContactAngleFvPatchScalarField::AngleCoefficient
(
    const scalar& x
) const
{
    return pow(pow(sin(x),3)/(2 - 3*cos(x) + pow(cos(x),3)),(1/3));
}


void Foam::saaDynamicKistlerAlphaContactAngleFvPatchScalarField::write(Ostream& os) const
{
    alphaContactAngleTwoPhaseFvPatchScalarField::write(os);
    os.writeKeyword("thetaA") << thetaA_ << token::END_STATEMENT << nl;
    os.writeKeyword("thetaR") << thetaR_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{
	makePatchTypeField
	(
		fvPatchScalarField,
		saaDynamicKistlerAlphaContactAngleFvPatchScalarField
	);

}

// ************************************************************************* //