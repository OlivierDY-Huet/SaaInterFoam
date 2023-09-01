Info << "Cell width calculation " << endl;

{
Info << "	Cell width" << endl;
vector meshDim_(runTime.controlDict().get<vector>("meshDim"));
scalar xCoef=meshDim_[0];
scalar yCoef=meshDim_[1];
scalar zCoef=meshDim_[2];

const faceList & ff = mesh.faces();
const pointField & pp = mesh.points();

label nMesh = 0;
scalar xDimSum = 0;
scalar yDimSum = 0;
scalar zDimSum = 0;

forAll(mesh.C(),celli)
{
	nMesh += 1;
    const cell & cc = mesh.cells()[celli];
    labelList pLabels(cc.labels(ff));
    pointField pLocal(pLabels.size(), vector::zero);
    forAll(pLabels, pointi)
	{
           pLocal[pointi] = pp[pLabels[pointi]];		   
	}
    scalar xDim = Foam::max(pLocal & vector(1,0,0)) - Foam::min(pLocal & vector(1,0,0));
	scalar yDim = Foam::max(pLocal & vector(0,1,0)) - Foam::min(pLocal & vector(0,1,0));
	scalar zDim = Foam::max(pLocal & vector(0,0,1)) - Foam::min(pLocal & vector(0,0,1));
	
	xDimSum += xDim;
	yDimSum += yDim;
	zDimSum += zDim;
	
	cellWidth[celli]= (xCoef*xDim + yCoef*yDim + zCoef*zDim)/(xCoef+yCoef+zCoef);
}
cellWidth.correctBoundaryConditions();

Info << "	Minimum cell width" << endl;
cellWidthMin.value() = gMin(cellWidth.field());


wallNeighb *= double(0);
const fvPatchList& patches = mesh.boundary();
forAll(patches, patchi)
{
	const fvPatch& curPatch = patches[patchi];
	if (isType<wallFvPatch>(curPatch))
	{	
		forAll(curPatch,facei)
		{
			label faceCelli = curPatch.faceCells()[facei];
			wallNeighb[faceCelli] = double(1);
		}			
	}
}
wallNeighb.correctBoundaryConditions();


}
