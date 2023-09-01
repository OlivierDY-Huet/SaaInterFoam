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

Application
    surfactantInterFoam.C

Group
    grpMultiphaseSolvers

Description
    Solver for two incompressible, isothermal immiscible fluids using a VOF
    (volume of fluid) phase-fraction based interface capturing approach,
    with optional mesh motion and mesh topology changes including adaptive
    re-meshing. One of the fluid containts surfactants.

\*---------------------------------------------------------------------------*/

#include <list>
#include <vector>
#include <string>
#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "fvcSmooth.H"

//-- Additional file --
#include "volPointInterpolation.H"
#include "interpolatePointToCell.H"	
#include "primitivePatchInterpolation.H"
//#include "pointMesh.H"
//#include "pointFields.H"
#include "fixedValuePointPatchFields.H"
#include "fvPatchFields.H"
#include "fvsPatchFields.H"
#include "alphaContactAngleTwoPhaseFvPatchScalarField.H"
#include "fixedValueFvsPatchFields.H"
#include "wallFvPatch.H" ///
//#include "addToRunTimeSelectionTable.H" //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Solver for two incompressible, isothermal immiscible fluids"
        " using VOF phase-fraction based interface capturing.\n"
        "With optional mesh motion and mesh topology changes including"
        " adaptive re-meshing."
    );
	
	bool iniMeshSwitch = false;
	bool iniSaaSwitch = false;

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
	
    #include "createFields.H"
	if (runTime.startTimeIndex()==0) //runTime.beginTime() runTime.startTime() runTime.startTimeIndex()
	{
		#include "setInitialAlpha.H"
		#include "setInitialLS.H"
		U *= double(0.0);
		mixture.correct();
		
		iniMeshSwitch = true;
		iniSaaSwitch  = false;
	}
	dimensionedScalar cellWidthMinOld = double(2)*cellWidthMin;
	//#include "cellWidth.H"
	//#include "meshLS.H"

    #include "createAlphaFluxes.H"
    #include "initCorrectPhi.H"
	#include "createUfIfPresent.H"


    turbulence->validate();

    if (!LTS)
    {
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"
    }
	
	#include "interfacePropertiesFields.H"
	
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;
    while (runTime.run())
    {
		

        #include "readDyMControls.H"
        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "CourantNo.H"
            #include "alphaCourantNo.H"
			#include "sigmaCourantNo.H"
			#include "sorptionCourantNo.H"
            #include "setDeltaT.H"
        }
		
		++runTime;
		
        Info<< "Time = " << runTime.timeName() << nl << endl;
		

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                mesh.update();
                if (mesh.changing())
				{
                    // Do not apply previous time-step mesh compression flux if the mesh topology changed
                    if (mesh.topoChanging())
                    {
                        talphaPhi1Corr0.clear();
                    }
                    gh = (g & mesh.C()) - ghRef;
                    ghf = (g & mesh.Cf()) - ghRef;
                    MRF.update();
                    if (correctPhi)
                    {
                        // Calculate absolute flux from the mapped surface velocity
                        phi = mesh.Sf() & Uf();
                        #include "correctPhi.H"
                        // Make the flux relative to the mesh motion
                        fvc::makeRelative(phi, U);
                        mixture.correct();
                    }
                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
					#include "cellWidth.H"
					#include "meshLS.H"
                }
								
            }
			
			if (iniMeshSwitch)
			{
				gh *= double(0.0);
				ghf *= double(0.0);
			}
			
			#include "alphaControls.H"
			#include "alphaEqnSubCycle.H"
			
            mixture.correct();		
			#include "Dirac.H"

            if (pimple.frozenFlow())
            {
                continue;
            }
			
			if (!iniMeshSwitch)
			{
				#include "saaConcEqn.H"
			}
			#include "saaConcInfArea.H"
			mixture.correct();
			
            #include "UEqn.H"
            // --- Pressure corrector loop
            while (pimple.correct())
            {	
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
			
        }
		
		////////////////
		
		Info<< "iniMeshSwitch: " << iniMeshSwitch  << endl;
		Info<< "iniSaaSwitch: " << iniSaaSwitch  << endl;
		
		
		if(iniSaaSwitch)
		{
			Info<< "	Add surfactant and velocity" << endl;
			#include "setInitialSaa.H"
			#include "setInitialU.H"
			
			mixture.correct();
			
			#include "initCorrectPhi.H"
			#include "createUfIfPresent.H"
			
			iniSaaSwitch  = false;
		}
		
		if(iniMeshSwitch)
		{
			Info<< "	Reset alpha" << endl;
			#include "setInitialAlpha.H"
			#include "setInitialLS.H"
			U *= double(0.0);		
			mixture.correct();
			#include "Dirac.H"
			
			#include "initCorrectPhi.H"
			#include "createUfIfPresent.H"
			
			Info<< "cellWidthMinOld: " << cellWidthMinOld  << endl;
			Info<< "cellWidthMin: " << cellWidthMin  << endl;
			
			if(cellWidthMin.value()<double(0.75)*cellWidthMinOld.value())
			{
				Info<< " Update cell width ref"  << endl;
				cellWidthMinOld = cellWidthMin;
			}
			else
			{
				Info<< " End of the alpha initialisation"  << endl;
				iniMeshSwitch = false;
				iniSaaSwitch = true;
			}					
		}
		
		////////////////
		
		
		#include "interfacePropertiesFields.H"
		
		runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
