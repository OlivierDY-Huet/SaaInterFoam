/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenFOAM Foundation
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

#include "saaConcDependentSurfaceTension.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceTensionModels
{
    defineTypeNameAndDebug(saaConcDependent, 0);
    addToRunTimeSelectionTable
    (
        surfaceTensionModel,
        saaConcDependent,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceTensionModels::saaConcDependent::saaConcDependent
(
    const dictionary& dict,
	const fvMesh& mesh
)
:
    surfaceTensionModel(mesh),
	cellWidthName_(dict.lookupOrDefault<word>("cellWidth", "cellWidth")),
    saaConcInfAreaName_(dict.lookupOrDefault<word>("saaConcInfArea", "saaConcInfArea"))//, 
	//sigma_(Function1<scalar>::New("sigma", dict))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surfaceTensionModels::saaConcDependent::~saaConcDependent()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::surfaceTensionModels::saaConcDependent::sigma() const
{

	// Initialisation
	Info << "	Initialisation" << endl;
	const volScalarField& cellWidth_ = mesh_.lookupObject<volScalarField>(cellWidthName_);
	const volScalarField& saaConcInfArea_ = mesh_.lookupObject<volScalarField>(saaConcInfAreaName_);
	
    auto tsigma = tmp<volScalarField>::New
    (
        IOobject
        (
            "sigma",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
		dimSigma,
		cellWidth_.boundaryField().types()
    );
    auto& sigma = tsigma.ref();
	
	// Surface tension calculation
	
	//	//Constants
	Info << "	Constants" << endl;
	dimensionedScalar cstR("cstR",dimensionSet(1,2,-2,-1,-1,0,0),scalar(8.314));
	dimensionedScalar cstT("cstT",dimensionSet(0,0,0,1,0,0,0),scalar(298.15));
	const dictionary& transportProperties = db().lookupObject<IOdictionary>("transportProperties");
	dimensionedScalar saaConcSurfMax("saaConcSurfMax",dimensionSet(0,-2,0,0,1,0,0),transportProperties);
	dimensionedScalar sigma0("sigma0",dimSigma,readScalar(transportProperties.subDict("sigma").lookup("sigma0")));
	dimensionedScalar saaN("saaN",dimensionSet(0,0,0,0,0,0,0),readScalar(transportProperties.subDict("sigma").lookup("saaN")));
	
	dimensionedScalar cstVal = saaN * cstR * cstT * saaConcSurfMax;
	
	Info << "	Calculation" << endl;
	sigma = sigma0 + cstVal * Foam::log (1- (saaConcInfArea_ / saaConcSurfMax));
	Info << "	End" << endl;
	/*
	scalar cstLim = min(saaConcSurfMax.value() , (1 - exp(-sigma0.value() / cstVal.value())) * saaConcSurfMax.value()) - VSMALL;
	
	
	// Correction of saaConcSurfMax if necessary
	Info << "	Correction (lim=" << cstLim << ")" << endl;
	volScalarField saaConcInfAreaCorr = saaConcInfArea_;
	
	
	forAll(saaConcInfArea_,celli)
	{
		if (saaConcInfArea_[celli] > cstLim)
		{
			Info << "Warning: the surfactant concentration is too high and will be corrected" << endl;
			saaConcInfAreaCorr[celli] = cstLim;
			const cell& c = mesh_.cells()[celli];
			forAll(c,facei)
			{
				label patchID = mesh_.boundaryMesh().whichPatch(c[facei]);
				if (patchID > -1)
				{
					label faceID = mesh_.boundaryMesh()[patchID].whichFace(c[facei]);
					if (saaConcInfAreaCorr.boundaryFieldRef()[patchID][faceID]> cstLim)
					{
						saaConcInfAreaCorr.boundaryFieldRef()[patchID][faceID]=cstLim;
					}
				}
			}
		}
	}
	*/
	/*
	int nCoorCell = 0;
	int nCoorPatch = 0;
	forAll(saaConcInfAreaCorr,celli)
	{
		if (saaConcInfAreaCorr[celli] > cstLim)
		{
			++nCoorCell;
			saaConcInfAreaCorr[celli] = cstLim;
		}
	}
	const volScalarField::Boundary& saaConcInfAreaCorrBf = saaConcInfAreaCorr.boundaryField();
	forAll(saaConcInfAreaCorrBf, patchi)
    {
		fvPatchScalarField& saaConcInfAreaCorrPatch = saaConcInfAreaCorr.boundaryFieldRef()[patchi];
		forAll(saaConcInfAreaCorrPatch, facei)
		if (saaConcInfAreaCorrPatch[facei]> cstLim)
		{
			++nCoorPatch;
			saaConcInfAreaCorrPatch[facei] = cstLim;
		}
    }
	
	if ((nCoorCell>0) || (nCoorPatch>0))
	{
		Info << "	Warning: the surfactant concentration. Nbr of corrections: cells=" << nCoorCell <<" ; patch faces=" << nCoorPatch <<endl;
	}
	
	Info << "	Calculation" << endl;
	
	sigma = sigma0 + cstVal * log (1- (saaConcInfAreaCorr / saaConcSurfMax));
	*/
	
	
    return tsigma;
}


bool Foam::surfaceTensionModels::saaConcDependent::readDict
(
    const dictionary& dict
)
{
    const dictionary& sigmaDict = surfaceTensionModel::sigmaDict(dict);
	
	cellWidthName_ = dict.lookupOrDefault<word>("cellWidth", "cellWidth");
	saaConcInfAreaName_ = dict.lookupOrDefault<word>("saaConcInfArea", "saaConcInfArea");
	//sigma_ = Function1<scalar>::New("sigma", sigmaDict);

    return true;
}


bool Foam::surfaceTensionModels::saaConcDependent::writeData
(
    Ostream& os
) const
{
    if (surfaceTensionModel::writeData(os))
    {
        //os  << sigma_() << token::END_STATEMENT << nl;
        return os.good();
    }

    return false;
}


// ************************************************************************* //
