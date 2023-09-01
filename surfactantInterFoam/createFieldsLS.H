 Info<<"\n LS variables" << endl;
 
bool LSIniSwitch = false;

//Psi

Info<<"   Creating field: " << "psi0"<< endl;
volScalarField psi0
(
	IOobject
	(
		"psi0",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ, //MUST_READ,
		IOobject::AUTO_WRITE//NO_WRITE
	),
	psi
);

/*
Info<<"   Creating field: " << "alpha1Sharp"<< endl;
volScalarField alpha1Sharp
(
	IOobject
	(
		"alpha1Sharp",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ, //MUST_READ,
		IOobject::AUTO_WRITE//NO_WRITE
	),
	alpha1
);
*/

Info<<"   Reading field: " << "psiInf"<< endl;
volScalarField psiInf
(
	IOobject
	(
		"psiInf",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	psi
);

Info<<"   Reading field: " << "psiInf"<< endl;
volScalarField psiRel
(
	IOobject
	(
		"psiRel",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	psi/cellWidthMin
);


Info<<"   Creating field: " << "innerInf"<< endl;
volScalarField innerInf
(
    IOobject
    (
        "innerInf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
	mesh,
	dimensionedScalar("dimScalar", dimensionSet(0,0,0,0,0,0,0), 0),
	alpha1.boundaryField().types()
);


Info<<"   Creating field: " << "signFunc"<< endl;
volScalarField signFunc
(
	IOobject
	(
		"signFunc",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE //NO_WRITE
	),
	mesh,
	dimensionedScalar("",dimless, -1),
	alpha1.boundaryField().types()
);



volVectorField gradPsi
(
	IOobject
	(
		"gradPsi",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	fvc::grad(psi)
);

volVectorField gradPsiWall
(
	IOobject
	(
		"gradPsiWall",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	gradPsi
);


///////////////////
/*
Info<<"   Creating field: " << "Gin"<< endl;
volScalarField Gin
(
	IOobject
	(
		"Gin",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
    mesh,
	dimensionedScalar("",dimless, 0.0),
	psi.boundaryField().types()
);

Info<<"   Creating field: " << "Gout"<< endl;
volScalarField Gout
(
	IOobject
	(
		"Gout",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
    mesh,
	dimensionedScalar("",dimless, 0.0),
	psi.boundaryField().types()
);

Info<<"   Creating field: " << "Gtot"<< endl;
volScalarField Gtot
(
	IOobject
	(
		"Gtot",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
    mesh,
	dimensionedScalar("",dimless, 0.0),
	psi.boundaryField().types()
);
*/

////////////
volScalarField aa
(
	IOobject
	(
		"aa",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ, //MUST_READ,
		IOobject::NO_WRITE//NO_WRITE
	),
	mag(fvc::grad(psi))
);

volScalarField bb
(
	IOobject
	(
		"bb",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ, //MUST_READ,
		IOobject::NO_WRITE//NO_WRITE
	),
	mag(fvc::grad(psi))
);

volScalarField cc
(
	IOobject
	(
		"cc",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ, //MUST_READ,
		IOobject::NO_WRITE//NO_WRITE
	),
	mag(fvc::grad(psi))
);

volScalarField dd
(
	IOobject
	(
		"dd",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ, //MUST_READ,
		IOobject::NO_WRITE//NO_WRITE
	),
	mag(fvc::grad(psi))
);

volScalarField ee
(
	IOobject
	(
		"ee",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ, //MUST_READ,
		IOobject::NO_WRITE//NO_WRITE
	),
	mag(fvc::grad(psi))
);

volScalarField ff
(
	IOobject
	(
		"ff",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ, //MUST_READ,
		IOobject::NO_WRITE//NO_WRITE
	),
	mag(fvc::grad(psi))
);

/////////////////////

volScalarField aaMod
(
	IOobject
	(
		"aaMod",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ, //MUST_READ,
		IOobject::NO_WRITE//NO_WRITE
	),
	mag(fvc::grad(psi))
);

volScalarField bbMod
(
	IOobject
	(
		"bbMod",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ, //MUST_READ,
		IOobject::NO_WRITE//NO_WRITE
	),
	mag(fvc::grad(psi))
);

volScalarField ccMod
(
	IOobject
	(
		"ccMod",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ, //MUST_READ,
		IOobject::NO_WRITE//NO_WRITE
	),
	mag(fvc::grad(psi))
);

volScalarField ddMod
(
	IOobject
	(
		"ddMod",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ, //MUST_READ,
		IOobject::NO_WRITE//NO_WRITE
	),
	mag(fvc::grad(psi))
);

volScalarField eeMod
(
	IOobject
	(
		"eeMod",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ, //MUST_READ,
		IOobject::NO_WRITE//NO_WRITE
	),
	mag(fvc::grad(psi))
);

volScalarField ffMod
(
	IOobject
	(
		"ffMod",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ, //MUST_READ,
		IOobject::NO_WRITE//NO_WRITE
	),
	mag(fvc::grad(psi))
);


/////////////////////


Info<<"   Reading field: " << "psiInf"<< endl;
volVectorField Ginv
(
	IOobject
	(
		"Ginv",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	fvc::grad(psi)
);


volVectorField Goutv
(
	IOobject
	(
		"Goutv",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	fvc::grad(psi)
);


volVectorField Gtotv
(
	IOobject
	(
		"Gtotv",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	fvc::grad(psi)
);
