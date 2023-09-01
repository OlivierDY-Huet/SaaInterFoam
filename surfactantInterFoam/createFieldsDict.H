Info<<"\n Read dictionaries for the saa and LS models" << endl;

//	//Transport properties
Info<<"   Reading transportProperties" << endl;
IOdictionary transportProperties
(
	IOobject
	(
		"transportProperties",
		runTime.constant(),
		mesh,
        IOobject::MUST_READ,//MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
	)
);

const dictionary& psiControls = mesh.solverDict(psi.name());

const dictionary& saaBulkControls = mesh.solverDict(saaConcBulk.name());
const dictionary& saaInfControls = mesh.solverDict(saaConcInf.name());
const dictionary& sigmaControls = mesh.solverDict(sigma.name());

