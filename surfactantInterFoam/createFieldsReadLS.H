Info<<"\n Read LS" << endl;

Info<<"   Reading field: " << "psi"<< endl;
volScalarField psi
(
	IOobject
	(
		"psi",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	mesh,
	dimensionedScalar("dimScalar", dimensionSet(0,1,0,0,0,0,0), 0),
	cellWidth.boundaryField().types()
);


// Coordonates

volVectorField XYZ_Coor
(
	IOobject
	(
		"XYZ_Coor",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	mesh,
	dimensionedVector("dimScalar", dimensionSet(0,1,0,0,0,0,0), vector(0.0,0.0,0.0)),
	cellWidth.boundaryField().types()
);
