Info<<"\n Read fields for the mesh" << endl;


Info << "Cell width " << endl;
volScalarField cellWidth
(
    IOobject
    (
        "cellWidth",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE //NO_WRITE
    ),
    mesh
);
dimensionedScalar cellWidthMin("cellWidthMin",cellWidth.dimensions(),scalar(0)); //dimensionedScalar cellWidthMin("cellWidthMin",dimensionSet(0,1,0,0,0,0,0),scalar(0));

/*
Info << "Cell width interpolated" << endl;
volScalarField cellWidthIntrpl
(
    IOobject
    (
        "cellWidthIntrpl",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE //NO_WRITE
    ),
    cellWidth
);
*/


Info<<"   Reading field: " << "wallNeighb"<< endl;
volScalarField wallNeighb
(
    IOobject
    (
        "wallNeighb",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE //NO_WRITE
    ),
	mesh,
	dimensionedScalar("",dimless, 0),
	cellWidth.boundaryField().types()
);
