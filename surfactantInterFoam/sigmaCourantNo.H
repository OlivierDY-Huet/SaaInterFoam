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

Global
    alphaCourantNo

Description
    Calculates and outputs the mean and maximum Courant Numbers.

\*---------------------------------------------------------------------------*/

scalar maxSigmaCo
(
    runTime.controlDict().get<scalar>("maxSigmaCo")
);

scalar sigmaCoNum = 0.0;

volScalarField saaSigma_ = mixture.sigmaIP();

if (mesh.nInternalFaces())
{
	volScalarField rho_ = alpha1*rho1 + (scalar(1.0)-alpha1)*rho2;
	
	scalarField tauSigmaRho
	(
		//pow(rho_.field()* mesh.V().field() / (saaSigma_.field() + VSMALL), 1/2)
		pow(rho_.field()* pow(cellWidth,3) / (saaSigma_.field() + VSMALL), 1/2)
	);
	
	scalarField tauSigmaMu
	(
		//mixture.mu() * pow(mesh.V().field(), 1/3) / (saaSigma_.field() + VSMALL)
		mixture.mu() * cellWidth / (saaSigma_.field() + VSMALL)
	);
		
	scalarField tauSigma = 0.5*
	(
	10 * tauSigmaMu 
	+ pow
		( 
		pow(10 * tauSigmaMu, 2)
		+ 4 * scalar(0.01) * pow(tauSigmaRho, 2)
		, 1/2
		)
	);
	sigmaCoNum  = runTime.deltaTValue()/(gMax(mixture.nearInterface()().primitiveField() * tauSigma)+SMALL);
}

Info<< "Additional Interface Courant Number max: " << sigmaCoNum << endl;

// ************************************************************************* //
