{
	
gradAlpha1 = fvc::grad(alpha1);	
gradPsi = fvc::grad(psi);
//volVectorField  gradAlpha1(fvc::grad(alpha1));

const volScalarField::Boundary& psibf = psi.boundaryField();
const fvPatchList& patches = mesh.boundary();
const polyBoundaryMesh& boundary = mesh.boundaryMesh();
forAll(boundary, patchi)
{
	const fvPatch& curPatch = patches[patchi];
	if (isA<fixedValueFvPatchScalarField>(psibf[patchi]))
	{
		
		/*		
		const fvPatchScalarField& alpha1p = alpha1.boundaryField()[patchi];
		scalarField alpha1pi = alpha1p.patchInternalField();
		scalarField condWall = pos0(alpha1pi - alpha1LimWall)*pos0(alpha2LimWall - alpha1pi);
		*/
		
		const fvPatchVectorField& gradAlpha1p = gradAlpha1.boundaryField()[patchi];
		vectorField gradAlpha1pi = gradAlpha1p.patchInternalField();
		vectorField nAlpha1pi = gradAlpha1pi/(mag(gradAlpha1pi)+SMALL);
		
		scalarField condWall = pos0(mag(gradAlpha1pi) * cellWidthMin.value() - relGradAlpha1LimWall); ///
			
		forAll(curPatch,facei)
		{
			label faceCelli = curPatch.faceCells()[facei];
			
			if(condWall[facei]>SMALL)
			{
				gradPsi[faceCelli] = nAlpha1pi[facei];
			}
			
		}
	}
}
gradPsi.correctBoundaryConditions();
}
