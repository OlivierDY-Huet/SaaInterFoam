{
	
Info << "Cell location " << endl;

/*
vector meshDim_(runTime.controlDict().get<vector>("meshDim"));
scalar xCoef=meshDim_[0];
scalar yCoef=meshDim_[1];
scalar zCoef=meshDim_[2];

forAll(mesh.cells(),celli)	
{
	X_Coor[celli]=mesh.C()[celli][0];
	Y_Coor[celli]=mesh.C()[celli][1];
	Z_Coor[celli]=mesh.C()[celli][2];
}

X_Coor.correctBoundaryConditions();
Y_Coor.correctBoundaryConditions();
Z_Coor.correctBoundaryConditions();
*/

XYZ_Coor.field() = mesh.C();
XYZ_Coor.correctBoundaryConditions();


Info << "epsilon and deltaTau " << endl;


epsilon = epsilonFact*cellWidthMin;
deltaTau = deltaTauFact*cellWidthMin;
}
