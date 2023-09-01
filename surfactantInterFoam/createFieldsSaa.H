 Info<<"\n Saa variables" << endl;


Info<<"   Creating field: " << "phic_"<< endl;
surfaceScalarField phic_
(
    IOobject
    (
        "phic_",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
	mag(phi/mesh.magSf())
);

Info<<"   Creating field: " << "phir_"<< endl;
surfaceScalarField phir_
(
    IOobject
    (
        "phir_",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
	mag(phi)
);


Info<<"   Creating field: " << "saaSink"<< endl;
volScalarField saaSink
(
    IOobject
    (
        "saaSink",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
	mesh,
	dimensionedScalar("dimScalar", dimensionSet(0,-3,-1,0,1,0,0), 0),
	cellWidth.boundaryField().types()
);

volScalarField saaSinkSp
(
    IOobject
    (
        "saaSinkSp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
	mesh,
	dimensionedScalar("dimScalar", dimensionSet(0,0,-1,0,0,0,0), 0),
	cellWidth.boundaryField().types()
);

volScalarField saaSinkSu
(
    IOobject
    (
        "saaSinkSu",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
	mesh,
	dimensionedScalar("dimScalar", dimensionSet(0,-3,-1,0,1,0,0), 0),
	cellWidth.boundaryField().types()
);

volScalarField saaSinkSpC
(
    IOobject
    (
        "saaSinkSpC",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
	mesh,
	dimensionedScalar("dimScalar", dimensionSet(0,-3,-1,0,1,0,0), 0),
	cellWidth.boundaryField().types()
);

Info<<"   Creating field: " << "saaSource"<< endl;
volScalarField saaSource
(
    IOobject
    (
        "saaSource",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
	mesh,
	dimensionedScalar("dimScalar", dimensionSet(0,-3,-1,0,1,0,0), 0),
	cellWidth.boundaryField().types()
);


volScalarField saaSourceSp
(
    IOobject
    (
        "saaSourceSp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
	mesh,
	dimensionedScalar("dimScalar", dimensionSet(0,0,-1,0,0,0,0), 0),
	cellWidth.boundaryField().types()
);

volScalarField saaSourceSu
(
    IOobject
    (
        "saaSourceSu",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
	mesh,
	dimensionedScalar("dimScalar", dimensionSet(0,-3,-1,0,1,0,0), 0),
	cellWidth.boundaryField().types()
);

volScalarField saaSourceSp1
(
    IOobject
    (
        "saaSourceSp1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
	mesh,
	dimensionedScalar("dimScalar", dimensionSet(0,0,-1,0,0,0,0), 0),
	cellWidth.boundaryField().types()
);

volScalarField saaSourceSp2
(
    IOobject
    (
        "saaSourceSp2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
	mesh,
	dimensionedScalar("dimScalar", dimensionSet(0,0,-1,0,0,0,0), 0),
	cellWidth.boundaryField().types()
);

volScalarField saaSourceSpC1
(
    IOobject
    (
        "saaSourceSpC1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
	mesh,
	dimensionedScalar("dimScalar", dimensionSet(0,-3,-1,0,1,0,0), 0),
	cellWidth.boundaryField().types()
);

volScalarField saaSourceSpC2
(
    IOobject
    (
        "saaSourceSpC2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
	mesh,
	dimensionedScalar("dimScalar", dimensionSet(0,-3,-1,0,1,0,0), 0),
	cellWidth.boundaryField().types()
);

Info<<"   Creating field: " << "wall sink"<< endl;
volScalarField saaWallSp
(
    IOobject
    (
        "saaWallSp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
	mesh,
	dimensionedScalar("dimScalar", dimensionSet(0,0,-1,0,0,0,0), 0),
	cellWidth.boundaryField().types()
);

volScalarField saaWallSpC
(
    IOobject
    (
        "saaWallSpC",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
	mesh,
	dimensionedScalar("dimScalar", dimensionSet(0,-3,-1,0,1,0,0), 0),
	cellWidth.boundaryField().types()
);


Info<<"   Creating field: " << "saaCBNew"<< endl;
volScalarField saaCBNew
(
    IOobject
    (
        "saaCBNew",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
	saaConcBulk
);

Info<<"   Creating field: " << "saaCINew"<< endl;
volScalarField saaCINew
(
    IOobject
    (
        "saaCINew",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
	saaConcInf
);

Info<<"   Creating field: " << "saaCWNew"<< endl;
volScalarField saaCWNew
(
    IOobject
    (
        "saaCWNew",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
	saaConcInf
);

/////////////////////////////////////////////////////
////////////////////////////////////////////////////




volVectorField cmpCB
(
    IOobject
    (
        "cmpCBr",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE//NO_WRITE
    ),
    mesh,
	dimensionedVector("dimScalar", dimensionSet(0,1,-1,0,0,0,0), vector (0,0,0)),
	cellWidth.boundaryField().types()
);

volVectorField cmpCI
(
    IOobject
    (
        "cmpCI",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE//NO_WRITE
    ),
    mesh,
	dimensionedVector("dimScalar", dimensionSet(0,1,-1,0,0,0,0), vector (0,0,0)),
	cellWidth.boundaryField().types()
);



///////////
Info<< "Field to write: saaConcInfIntrpl" << endl;
volScalarField saaConcInfIntrpl
(
    IOobject
    (
        "saaConcInfIntrpl",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE//NO_WRITE
    ),
    saaConcInf
);


Info<< "Field to write: gradAlpha1" << endl;
volVectorField gradAlpha1
(
    IOobject
    (
        "gradAlpha1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE//NO_WRITE
    ),
    fvc::grad(alpha1)
);


Info<< "Field to write: Dirac" << endl;
volScalarField DiracInf
(
    IOobject
    (
        "DiracInf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE//NO_WRITE
    ),
    mesh,
	dimensionedScalar("dimScalar",dimensionSet(0,-1,0,0,0,0,0), 0),
	cellWidth.boundaryField().types()
);

Info<< "Field to write: DiracWall" << endl;
volScalarField DiracWall
(
    IOobject
    (
        "DiracWall",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE//NO_WRITE
    ),
    mesh,
	dimensionedScalar("dimScalar",dimensionSet(0,-1,0,0,0,0,0), 0),
	cellWidth.boundaryField().types()
);

Info<< "Field to write: DiracIntrpl" << endl;
volScalarField DiracIntrpl
(
    IOobject
    (
        "DiracIntrpl",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE//NO_WRITE
    ),
    mesh,
	dimensionedScalar("dimScalar",DiracInf.dimensions(), 0),
	cellWidth.boundaryField().types()
);

Info<< "Field to write: coefMagGrad" << endl;
volScalarField coefMagGrad
(
    IOobject
    (
        "coefMagGrad",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE//NO_WRITE
    ),
    mesh,
	dimensionedScalar("dimScalar",dimensionSet(0,0,0,0,0,0,0), 0),
	cellWidth.boundaryField().types()
);

Info<< "Field to write: Heaviside" << endl;
volScalarField HeavisideInf
(
    IOobject
    (
        "HeavisideInf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE//NO_WRITE
    ),
    mesh,
	dimensionedScalar("dimScalar",dimensionSet(0,0,0,0,0,0,0), 0),
	cellWidth.boundaryField().types()
);



Info<<"   Creating field: " << "cellDistInf"<< endl;
scalar maxCellDist = 6;
volScalarField cellDistInf
(
    IOobject
    (
        "cellDistInf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE//NO_WRITE
    ),
	mesh,
	dimensionedScalar("dimScalar", dimensionSet(0,0,0,0,0,0,0), 0),
	cellWidth.boundaryField().types()
);

Info<< "Field to check: checkInfLink" << endl;
volScalarField checkInfLink
(
    IOobject
    (
        "checkInfLink",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE//NO_WRITE
    ),
	mesh,
	dimensionedScalar("dimScalar", dimensionSet(0,0,0,0,0,0,0), 0),
	cellWidth.boundaryField().types()
);



volScalarField saaConcSurfRel
(
    IOobject
    (
        "saaConcSurfRel",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE//NO_WRITE
    ),
	mesh,
	dimensionedScalar("dimScalar", dimensionSet(0,0,0,0,0,0,0), 0),
	cellWidth.boundaryField().types()
);


////////////////


volVectorField coefInfGrad
(
    IOobject
    (
        "coefInfGrad",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE//NO_WRITE
    ),
	mesh,
	dimensionedVector("dimScalar", dimensionSet(0,0,0,0,0,0,0), vector (0,0,0)),
	cellWidth.boundaryField().types()
);

volVectorField coefBulkGrad
(
    IOobject
    (
        "coefBulkGrad",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE//NO_WRITE
    ),
	mesh,
	dimensionedVector("dimScalar", dimensionSet(0,0,0,0,0,0,0), vector (0,0,0)),
	cellWidth.boundaryField().types()
);

volVectorField coefInfLim
(
    IOobject
    (
        "coefInfLim",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE//NO_WRITE
    ),
	mesh,
	dimensionedVector("dimScalar", dimensionSet(0,0,0,0,0,0,0), vector (0,0,0)),
	cellWidth.boundaryField().types()
);

volVectorField coefBulkLim
(
    IOobject
    (
        "coefBulkLim",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE//NO_WRITE
    ),
	mesh,
	dimensionedVector("dimScalar", dimensionSet(0,0,0,0,0,0,0), vector (0,0,0)),
	cellWidth.boundaryField().types()
);

////////

Info<<"   Creating field: " << "saaConcInfAreaSharp"<< endl;
volScalarField saaConcInfAreaSharp
(
    IOobject
    (
        "saaConcInfAreaSharp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
	saaConcInfArea
);


volScalarField sorptBarrierInf
(
    IOobject
    (
        "sorptBarrierInf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE//NO_WRITE
    ),
	mesh,
	dimensionedScalar("dimScalar", dimensionSet(0,0,0,0,0,0,0), 0),
	cellWidth.boundaryField().types()
);


