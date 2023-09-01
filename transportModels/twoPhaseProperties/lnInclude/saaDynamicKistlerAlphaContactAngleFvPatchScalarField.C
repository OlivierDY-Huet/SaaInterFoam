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
	
	/*
	const dictionary& transportProperties = db().lookupObject<IOdictionary>("transportProperties");
	scalar nu1(readScalar(transportProperties.subDict("water").lookup("nu")));
	scalar rho1(readScalar(transportProperties.subDict("water").lookup("rho")));
	scalar mup = nu1* rho1;
	*/
	//Info << db().names() << endl;
	
	//const fvPatchField<scalar>& mup = patch().lookupPatchField<volScalarField, scalar>("mu");
	//const fvPatchField<scalar>& sigmap = patch().lookupPatchField<volScalarField, scalar>("sigmaIP");
	//const fvPatchField<vector>& UpSmoothp = patch().lookupPatchField<volVectorField, vector>("UKistler");
	//const fvPatchField<scalar>& alphap = patch().lookupPatchField<volScalarField, scalar>("alpha.water");
	//const fvsPatchField<vector>& nHatKistlerp = patch().lookupPatchField<surfaceVectorField, vector>("nHatKistler");
	
	//const polyPatch& pp = alphap.patch().patch();
	
	const fvPatchField<scalar>& sigmap = patch().lookupPatchField<volScalarField, scalar>("sigmaIP");
	const fvPatchField<scalar>& mup = patch().lookupPatchField<volScalarField, scalar>("muEff");
	
	
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

    //eb - Calculate InverseHoffmanFunction for thetaA and thetaR using
    // RiddersRoot
    RiddersRoot RRInvHoffFuncThetaA(InvHoffFuncThetaA, 1.e-10);
    scalar InvHoffFuncThetaAroot = RRInvHoffFuncThetaA.root(0,65);

    RiddersRoot RRInvHoffFuncThetaR(InvHoffFuncThetaR, 1.e-10);
    scalar InvHoffFuncThetaRroot = RRInvHoffFuncThetaR.root(0,65);

    //eb - Calculate and return the value of contact angle on patch faces,
    //     a general approach: the product of Uwall and nWall is negative
    //     for advancing and positiv for receding motion.
    //     thetaDp is initialized to the equilibrium contact angle https://pubs.acs.org/doi/10.1021/la049410h
	
	scalar coefA = AngleCoefficient(convertToRad*thetaA_);
	scalar coefR = AngleCoefficient(convertToRad*thetaR_);
	scalar theta0 = constant::mathematical::pi/2; //= acos((coefA*cos(convertToRad*thetaA_) + coefR*cos(convertToRad*thetaR_))/(coefA+coefR));
	Info << "Equilibrium CA: " << convertToDeg*theta0 << " [deg]" << endl;
	
    scalarField thetaDp(patch().size(), convertToRad*theta0);
    forAll(uwall, pfacei)
    {
		//if (alphap[pfacei] < scalar(0.99) && alphap[pfacei] > scalar(0.01))
		//{
			scalar thetaK = convertToRad*theta0;
			if(uwall[pfacei] < 0.0)
			{
				thetaK = HoffmanFunction(   Ca[pfacei]
												   + InvHoffFuncThetaAroot);
			}
			else if (uwall[pfacei] > 0.0)
			{
				thetaK = HoffmanFunction(   max(Ca[pfacei],-InvHoffFuncThetaRroot+SMALL)
												   + InvHoffFuncThetaRroot);
			}
			
			thetaDp[pfacei] = min(max(thetaK,scalar(0.0)+convertToRad*SMALL),constant::mathematical::pi-convertToRad*SMALL);
		//}
		//else
		//{
		//	thetaDp[pfacei] = constant::mathematical::pi/2;
		//}
		//Info << "	Patch corr: alpha = " << alphap[pfacei] << "[-] ; ST = " << sigmap[pfacei] << " [N/m] ; Uwall = " << uwall[pfacei] <<" [m/s] ; mu = " << mup[pfacei] << "[kg/m s] Ca = " << Ca[pfacei] << " [-] ; CA = " << convertToDeg*thetaDp[pfacei] << " [deg]"<< endl;
    	//Info << "	Patch corr: alpha = " << alphap[pfacei] << "[-] ; ST = " << sigmap[pfacei] << " [N/m] ; Uwall = " << uwall[pfacei] <<" [m/s] ; mu = " << mup << "[kg/m s] Ca = " << Ca[pfacei] << " [-] ; CA = " << convertToDeg*thetaDp[pfacei] << " [deg]"<< endl;
		//Info << "	Patch corr: alpha = " << alphap[pfacei] << "[-] ; ST = " << sigmap[pfacei] << " [N/m] ; Uwall = " << Uwall[pfacei] <<" [m/s] ; nWall = " << nWall[pfacei] <<" [-] ; uwall = " << uwall[pfacei] <<" [m/s] ; mu = " << mup << "[kg/m s] Ca = " << Ca[pfacei] << " [-] ; CA = " << convertToDeg*thetaDp[pfacei] << " [deg]"<< endl;
	}
	
	// Correction
	/*
	forAll(uwall, pfacei)
	{
		//std::vector<std::list<int> > faceDistInf;
		
		std::vector<int> nFaces;

		const labelList& faceEdges = pp.faceEdges()[pfacei];
		forAll(faceEdges, ppedgei)
		{
			const labelList& neigbourFaces = pp.edgeFaces()[faceEdges[ppedgei]];
			forAll(neigbourFaces,ppfacei)
			{
				nFaces.push_back(neigbourFaces[ppfacei]);
			}
		}	
		int ixdInf = 0;
		forAll(nFaces,ppfacei)
		{
			if (alphap[nFaces[ppfacei]] >0.5)
			{
				if (alphap[nFaces[ppfacei]]<alphap[nFaces[ixdInf]])
				{
					ixdInf = ppfacei;
				}
			}
		}
		
		if (( alphap[nFaces[ixdInf]] > scalar(0.01) ) && ( alphap[nFaces[ixdInf]] < scalar(0.99) ))
		{
			thetaDp[pfacei] = thetaDp[nFaces[ixdInf]];
		}
		Info << "	Patch corr: alpha = " << alphap[pfacei] << "[-] ; ST = " << sigmap[pfacei] << " [N/m] ; Uwall = " << uwall[pfacei] <<" [m/s] ; mu = " << mup[pfacei] << "[kg/m s] Ca = " << Ca[pfacei] << " [-] ; CA = " << convertToDeg*thetaDp[pfacei] << " [deg]"<< endl;
		
		
	}
	*/
    
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