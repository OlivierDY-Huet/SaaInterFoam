{

volVectorField upwindDir(fvc::reconstruct(SSF));

const volScalarField::Boundary& psiDummybf = psiDummy.boundaryField();
const volScalarField::Boundary& psibf = psi.boundaryField();
const fvPatchList& patches = mesh.boundary();
const polyBoundaryMesh& boundary = mesh.boundaryMesh();
forAll(boundary, patchi)
{
	if (isA<fixedValueFvPatchScalarField>(psibf[patchi]))
	{
		
		fixedValueFvPatchScalarField& psip = 
		const_cast<fixedValueFvPatchScalarField&>
        (
            refCast<const fixedValueFvPatchScalarField>
            (
                psibf[patchi]
            )
        );
		
		scalarField psipi = psip.patchInternalField();
		
		fixedValueFvPatchScalarField& psiDummyp = 
		const_cast<fixedValueFvPatchScalarField&>
        (
            refCast<const fixedValueFvPatchScalarField>
            (
                psiDummybf[patchi]
            )
        );
		
		
		const fvPatchVectorField& upwindDirp = upwindDir.boundaryField()[patchi];
		vectorField upwindDirpi = upwindDirp.patchInternalField();
		vectorField nUW = upwindDirpi/(mag(upwindDirpi)+SMALL);
				
		const fvPatch& curPatch = patches[patchi];
		vectorField nf = curPatch.nf();
		
		scalarField coef = pos0(nUW & -nf);
		
		psiDummyp == coef * psip + (1.0-coef) * psipi;
		psiDummyp.evaluate();
	}
}

psiDummy.correctBoundaryConditions();

gradPsiDummy = fvc::grad(psiDummy);
}
