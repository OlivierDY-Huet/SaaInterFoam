Info<< "Field to check: sigma" << endl;

Info<<"   Reading transportProperties: " << "saaConcSurfEq" << endl;
dimensionedScalar saaConcSurfEq
(
	"saaConcSurfEq",
	dimensionSet(0,-2,0,0,1,0,0),
	saaConcSurfMax.value()/2
);

Info<<"   Reading transportProperties: " << "saaConcSurfEqMax" << endl;
dimensionedScalar saaConcSurfEqMax
(
	"saaConcSurfEqMax",
	dimensionSet(0,-2,0,0,1,0,0),
	saaConcSurfMax.value()/2
);

Info<<"   Reading transportProperties: " << "saaConcSurfIni" << endl;
dimensionedScalar saaConcSurfIni
(
	"saaConcSurfIni",
	dimensionSet(0,-2,0,0,1,0,0),
	transportProperties
);


Info<< "Field to check: K" << endl;
volScalarField K
(
    IOobject
    (
        "K",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE//NO_WRITE
    ),
	mesh,
	dimensionedScalar("dimScalar", dimensionSet(0,-1,0,0,0,0,0), 0),
	cellWidth.boundaryField().types()
);


Info<< "Field to check: gradSigma" << endl;
volVectorField gradSigma
(
    IOobject
    (
        "gradSigma",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
	fvc::grad(sigma)
);


Info<< "Field to check: gradHeaviside" << endl;
volVectorField gradHeaviside
(
    IOobject
    (
        "gradHeaviside",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
	fvc::grad(HeavisideInf)
);

Info<< "Field to check: nHat" << endl;
volVectorField nHatInf
(
    IOobject
    (
        "nHatInf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::grad(alpha1)/(mag(fvc::grad(alpha1))+dimensionedScalar("",dimensionSet(0,-1,0,0,0,0,0),1)*SMALL)//fvc::reconstruct(mixture.nHatf()* mesh.magSf())
);


Info<< "Field to check: SurfaceTensionForce" << endl;
volVectorField SurfaceTensionForce
(
    IOobject
    (
        "SurfaceTensionForce",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::grad(alpha1)*mag(fvc::grad(sigma))
);

Info<< "Field to check: MarangoniForce" << endl;
volVectorField MarangoniForce
(
    IOobject
    (
        "MarangoniForce",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mag(fvc::grad(alpha1))*fvc::grad(sigma)
);

Info<< "Field to check: scaleCoefInf" << endl;
volScalarField scaleCoefInf
(
    IOobject
    (
        "scaleCoefInf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    alpha1,
	cellWidth.boundaryField().types()
);



Info<< "Field to check: MarangoniForceMod" << endl;
volVectorField MarangoniScaledForce
(
    IOobject
    (
        "MarangoniScaledForce",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mag(fvc::grad(alpha1))*fvc::grad(sigma)
);


Info<< "Field to check: DiracInfWall" << endl;
volScalarField DiracInfWall
(
    IOobject
    (
        "DiracInfWall",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mag(fvc::grad(alpha1))
);

//////////////////////////////////////////////////////////////


Info<<"   Reading transportProperties: " << "phisSwitch" << endl;
const int phisSwitch = transportProperties.lookupOrDefault<int>("phisSwitch", 0);


volVectorField phirField
(
    IOobject
    (
        "phirField",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    U
);



Info<< "Field to check: phisBulkField" << endl;
volVectorField phisBulkField
(
    IOobject
    (
        "phisBulkField",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    U
);


