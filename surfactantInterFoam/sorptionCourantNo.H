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

scalar maxSorptionCo
(
    runTime.controlDict().get<scalar>("maxSorptionCo")
);

scalar sorptionCoNum = 0.0;

if (coefAd.value() > 0.0)
{
if (mesh.nInternalFaces())
{
	
	scalarField highSorptionCell(pos0(mixture.Dirac().field()*epsilon.value()-double(0.25)));
	//sorptionCoNum = gMax(mixture.nearInterface2Epsilon()().primitiveField() * mag(saaSource.field()) * cellWidth.field() /saaConcSurfEq.value()) * runTime.deltaTValue();
	//sorptionCoNum = gMax(mixture.nearInterface2Epsilon()().primitiveField() * mag(saaSource.field()) / (Dirac.field() * saaConcSurfEq.value() + VSMALL)) * runTime.deltaTValue();
	//sorptionCoNum = gMax(mixture.nearInterface1Epsilon()().primitiveField() * mag(saaSource.field()) / (mixture.Dirac().field() * saaConcSurfEq.value() + VSMALL)) * runTime.deltaTValue();
	sorptionCoNum = gMax(highSorptionCell * mag(saaSource.field()) / ((mixture.Dirac().field() + coefInfCorrSMALL/epsilon.value()) * saaConcSurfMax.value())) * runTime.deltaTValue();
}
}

Info<< "Additional Sorption Courant Number max: " << sorptionCoNum << endl;

// ************************************************************************* //
