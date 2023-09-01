Info<< "   Dirac " <<endl;

// Heaviside
HeavisideInf = mixture.Heaviside();
HeavisideInf.correctBoundaryConditions();

// Dirac inf
DiracInf = mixture.Dirac();
DiracInf.correctBoundaryConditions();

innerInf = mixture.nearInterface1Epsilon(); //pos0(DiracInf*epsilon-SMALL);
innerInf.correctBoundaryConditions();


// Dirac wall
if (saaWallEffCoef>SMALL)
{
	DiracWall *= double(0);
	const fvPatchList& patches = mesh.boundary();
	forAll(patches, patchi)
	{
		const fvPatch& curPatch = patches[patchi];
		if (isType<wallFvPatch>(curPatch))
		{	
			forAll(curPatch,facei)
			{
				label faceCelli = curPatch.faceCells()[facei];
				DiracWall[faceCelli] = mag(scaleCoefInf[faceCelli])/(double(2)*cellWidth[faceCelli]);
			}			
		}
	}
}

