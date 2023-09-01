#include "createFieldsMesh.H"
#include "cellWidth.H"

#include "createFieldsReadLS.H"
#include "createFieldsReadSaa.H"

#include "createFieldsDict.H"


#include "createFieldsOptionsDrop.H"
#include "createFieldsOptionsLS.H"
#include "createFieldsOptionsSaa.H"

#include "createFieldsCA.H"


Info<<"\n Interfoam variables" << endl;

#include "createRDeltaT.H"

Info<<"   Reading field: " << "p_rgh"<< endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<<"   Reading field: " << "U"<< endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"

Info<<"   Phase properties"<< endl;
immiscibleIncompressibleTwoPhaseMixture mixture(U, phi, psi, epsilon);

Info<<"   Phase properties: " << "alpha 1 and 2"<< endl;
volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

Info<<"   Phase properties: " << "rho 1 and 2"<< endl;
const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();

Info<<"   Phase properties: " << "nu 1 and 2"<< endl;
const dimensionedScalar& nu1 = mixture.nu1(); 
const dimensionedScalar& nu2 = mixture.nu2();

Info<<"   Creating field: " << "rho"<< endl;
volScalarField rho // Need to store rho for ddt(rho, U)
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    alpha1*rho1 + alpha2*rho2
);
rho.oldTime();

Info<<"   Creating field: " << "rhoPhi"<< endl;
surfaceScalarField rhoPhi // Mass flux
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);


// Construct incompressible turbulence model
Info<<"   Turbulent model"<< endl;
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, mixture)
);

Info<<"   Gravity"<< endl;
#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"

Info<<"   Creating field: " << "p"<< endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh;
}

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());

// MULES compressed flux is registered in case scalarTransport FO needs it.
Info<<"   Creating field: " << "alphaPhiUn"<< endl;
surfaceScalarField alphaPhiUn
(
    IOobject
    (
        "alphaPhiUn",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(phi.dimensions(), Zero)
);

#include "createMRF.H"
#include "createFvOptions.H"

#include "createFieldsSaa.H"
#include "createFieldsLS.H"
#include "createFieldsIP.H"


Info<<"   Creating field: " << "END"<< endl;
