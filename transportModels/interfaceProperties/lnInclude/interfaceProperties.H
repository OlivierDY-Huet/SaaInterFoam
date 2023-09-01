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

Class
    Foam::interfaceProperties

Description
    Contains the interface properties.

    Properties to aid interFoam:
    -# Correct the alpha boundary condition for dynamic contact angle.
    -# Calculate interface curvature.

SourceFiles
    interfaceProperties.C

\*---------------------------------------------------------------------------*/

#ifndef interfaceProperties_H
#define interfaceProperties_H

#include "IOdictionary.H"
#include "surfaceTensionModel.H"
#include "volFields.H"
#include "surfaceFields.H"

//#include "volPointInterpolation.H"
//#include "interpolatePointToCell.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class interfaceProperties Declaration
\*---------------------------------------------------------------------------*/

class interfaceProperties
{
    // Private data

        //- Keep a reference to the transportProperties dictionary
        const dictionary& transportPropertiesDict_;

        //- Compression coefficient
        scalar cAlpha_;

        //- Surface tension
        autoPtr<surfaceTensionModel> sigmaPtr_;
		
		//- Input fields
		const volScalarField& alpha1_;
        const volVectorField& U_;
		
		//- LS inputs
		const volScalarField& psi_;
		const dimensionedScalar& epsilon_;
		
		//- Heaviside
		volScalarField Heaviside_;
		surfaceScalarField Heavisidef_;
		
		//- Dirac
		volScalarField Dirac_;
		surfaceScalarField Diracf_;
		
        //- Curvature 
		const dimensionedScalar deltaN_;
		
		volVectorField nHat_;
		surfaceVectorField nHatAlphafv_;
		surfaceVectorField nHatHeavisidefv_;
		surfaceVectorField nHatfv_;
		
		surfaceScalarField compStrCoeff_;
		
		volScalarField K_;
		
		//- Contact angle
		volScalarField thetaRad_;
		
		//- Surface tension
		volScalarField sigmaIP_;
		
		
    // Private Member Functions

        //- No copy construct
        interfaceProperties(const interfaceProperties&) = delete;

        //- No copy assignment
        void operator=(const interfaceProperties&) = delete;

        //- Correction for the boundary condition on the unit normal nHat on
        //  walls to produce the correct contact dynamic angle
        //  calculated from the component of U parallel to the wall
        void correctContactAngle
        (
            surfaceVectorField::Boundary& nHat,
            const surfaceVectorField::Boundary& gradAlphaf,
			volScalarField& thetaRad_
        ) const;
		
		
		void smoothen
		(
			volScalarField& smooth_func
		) const;
		
		
		//- Re-calculate Heaviside
		void calculateHeaviside();
		
		//- Re-calculate Dirac
		void calculateDirac();

        //- Re-calculate the interface curvature
        void calculateK();
		
		//- Re-calculate ST
		void calculateSurfaceTension();

public:

    // Constructors

        //- Construct from volume fraction field gamma and IOdictionary
        interfaceProperties
        (
            const volScalarField& alpha1,
            const volVectorField& U,
			const volScalarField& psi,
			const dimensionedScalar& epsilon,
            const IOdictionary&
        );


    // Member Functions

        scalar cAlpha() const
        {
            return cAlpha_;
        }
		
        const dimensionedScalar& deltaN() const
        {
            return deltaN_;
        }
		
		volVectorField nHat() const
        {
            return nHat_;
        }
		
		volScalarField Heaviside() const
        {
            return Heaviside_;
        }
		
		volScalarField Dirac() const
        {
            return Dirac_;
        }	
		
		surfaceVectorField nHatAlphafv() const
        {
            return nHatAlphafv_;
        }
		
		surfaceVectorField nHatHeavisidefv() const
        {
            return nHatHeavisidefv_;
        }
		
		surfaceVectorField nHatfv() const
        {
            return nHatfv_;
        }
		
		surfaceScalarField compStrCoeff() const
        {
            return compStrCoeff_;
        }
		
		volScalarField KCoef() const
        {
            return K_;
        }
		
		volScalarField thetaRad() const
        {
            return thetaRad_;
        }
			
		volScalarField sigmaIP() const
        {
            return sigmaIP_;
        }
		
		
        tmp<volScalarField> sigmaK() const;

        tmp<volVectorField> surfaceTensionForce() const;
		
		tmp<volVectorField> surfaceTensionDensityScaledForce() const;
		
		tmp<surfaceScalarField> surfaceTensionForcef() const;
		
		tmp<surfaceScalarField> surfaceTensionDensityScaledForcef() const;
		
		tmp<volVectorField> MarangoniForce() const;
		
		tmp<volVectorField> MarangoniDensityScaledForce() const;
		
		tmp<surfaceScalarField> MarangoniForcef() const;
		
		tmp<surfaceScalarField> MarangoniDensityScaledForcef() const;
		
        //- Indicator of the proximity of the interface
        //  Field values are 1 near and 0 away for the interface.
        tmp<volScalarField> nearInterface() const;
		
		tmp<volScalarField> nearInterface1Epsilon() const;
		
        void correct();

        //- Read transportProperties dictionary
        bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
