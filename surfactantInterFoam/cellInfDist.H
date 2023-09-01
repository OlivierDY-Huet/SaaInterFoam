Info <<"Cell distance field creation"<<endl;

Info <<"	Shared variables"<<endl;
std::vector<std::list<int> > cellDistInfList;
std::list<int> cellDistInfLvl;
label meshLvl = int(0);
scalar lim1lpha = double(0.5);

Info <<"	Initialisation"<<endl;
const labelListList& cellCells = mesh.cellCells();
std::list<int> cellDistInfLvlPrev;


Info <<"	Interface"<<endl;

volScalarField nearInterface = innerInf;//mixture.nearInterface1Epsilon();
cellDistInf = double(1)-nearInterface;
forAll(cellDistInf,celli)
{
	if(nearInterface[celli]>0)
	{
		cellDistInf[celli]=double(meshLvl);
		cellDistInfLvl.push_back(celli);
	}
	else
	{
		cellDistInf[celli]=sign(alpha1[celli]-lim1lpha)*maxCellDist;
	}
}

cellDistInf.correctBoundaryConditions();

cellDistInfLvl.unique();
cellDistInfList.push_back(cellDistInfLvl);
cellDistInfLvl.clear();


Info <<"   Neighbours"<<endl;
//	//Lvl >0

for (meshLvl=1; meshLvl<maxCellDist; ++meshLvl)
{
	//Evaluation
	cellDistInfLvlPrev=cellDistInfList[meshLvl-1];
	for (const int& celli : cellDistInfLvlPrev)
	{
		const labelList& nNeighbourCells = cellCells[celli];
		forAll(nNeighbourCells, ni)
		{	
			int cellNi=nNeighbourCells[ni];
			if (abs(cellDistInf[cellNi])>(meshLvl))
			{					
				cellDistInfLvl.push_back(cellNi);
				cellDistInf[cellNi] = sign(alpha1[cellNi]-lim1lpha)*meshLvl;
			}
		}
	}
	//Correction processor 
	//cellDistInf.correctBoundaryConditions();
	const fvPatchList& patches = mesh.boundary();
	forAll(patches, patchi)
	{
		const fvPatch& curPatch = patches[patchi];
		if (isType<processorFvPatch>(curPatch))
		{
			const fvPatchScalarField& cellDistInfBf = cellDistInf.boundaryField()[patchi];
			forAll(curPatch,facei)
			{	
				int celli = curPatch.faceCells()[facei];			
				if ((mag(cellDistInfBf[facei])<(meshLvl)) && (mag(cellDistInf[celli])>meshLvl))
				{	
					cellDistInf[celli]=sign(alpha1[celli]-lim1lpha)*meshLvl;
					cellDistInfLvl.push_back(celli);
				}
			}
		}
	}
	cellDistInf.correctBoundaryConditions();
	
	// List transfer
	cellDistInfLvl.unique();
	cellDistInfList.push_back(cellDistInfLvl);
	cellDistInfLvl.clear();
}

