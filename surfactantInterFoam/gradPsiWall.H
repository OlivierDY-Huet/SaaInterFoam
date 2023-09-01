{
	
//reset references;
gradPsiWall *= double(0);

//Modify references at wall;
const volScalarField::Boundary& abf = alpha1.boundaryField();
const fvPatchList& patches = mesh.boundary();
const polyBoundaryMesh& boundary = mesh.boundaryMesh();
forAll(boundary, patchi)
{
	if (isA<alphaContactAngleTwoPhaseFvPatchScalarField>(abf[patchi]))
	{
		
		const fvPatch& curPatch = patches[patchi];
		vectorField nf = curPatch.nf();
		
		const fvPatchVectorField& gradPsip = gradPsi.boundaryField()[patchi];
		vectorField gradPsipi = gradPsip.patchInternalField();
		vectorField nPsi = gradPsipi/(mag(gradPsipi)+SMALL);
		
		vectorField gradPsiWallpi = gradPsipi * mag(nf & nPsi);
		
		forAll(curPatch,facei)
		{	
			label faceCelli = curPatch.faceCells()[facei];	

			gradPsiWall[faceCelli] = gradPsiWallpi[facei];
		}
		
		
	}
}
gradPsiWall.correctBoundaryConditions();

}
