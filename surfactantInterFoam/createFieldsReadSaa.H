Info<<"\n Read Saa" << endl;


//////////////////



//	//Fields

Info<<"   Reading field: " << "saaConcBulk"<< endl;
volScalarField saaConcBulk
(
    IOobject
    (
        "saaConcBulk",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<<"   Reading field: " << "saaConcInf"<< endl;
volScalarField saaConcInf
(
    IOobject
    (
        "saaConcInf",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<<"   Reading field: " << "saaConcWall"<< endl;
volScalarField saaConcWall
(
    IOobject
    (
        "saaConcWall",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<<"   Creating field: " << "saaConcInfArea"<< endl;
volScalarField saaConcInfArea
(
    IOobject
    (
        "saaConcInfArea",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
	mesh,
	dimensionedScalar("",dimensionSet(0,-2,0,0,1,0,0), 0),
	cellWidth.boundaryField().types()
);


volScalarField sigma
(
    IOobject
    (
        "sigma",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
	mesh,
	dimensionedScalar("",dimensionSet(1,0,-2,0,0,0,0), 0),
	cellWidth.boundaryField().types()
);
