Info<<"   Reading field: " << "thetaDeg"<< endl;
volScalarField thetaDeg
(
    IOobject
    (
        "thetaDeg",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
	dimensionedScalar("",dimless, 0),
	cellWidth.boundaryField().types()
);
