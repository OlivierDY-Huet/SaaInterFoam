Info<< "\nArea specific surface concentration" << endl;
{


Info << "   Interpolation" << endl;
const dimensionedScalar dt("dt",dimensionSet(0,0,1,0,0,0,0),mesh.time().deltaTValue());

coefMagGrad = innerInf * mag(fvc::grad(alpha1))/(mag(fvc::grad(HeavisideInf))+SMALL/cellWidthMin);

saaConcInfIntrpl = saaConcInf * innerInf;
saaConcInfIntrpl.correctBoundaryConditions();

DiracIntrpl = DiracInf * innerInf;
DiracIntrpl.correctBoundaryConditions();

//saaConcInfAreaSharp = saaConcInfIntrpl / (DiracIntrpl + SMALL*dimensionedScalar("",dimensionSet(0,-1,0,0,0,0,0),1.0));

Info << "   -" << endl;

saaConcInfAreaSharp = innerInf * saaConcInf / (DiracInf+coefInfCorrSMALL/epsilon);

Info << "   -" << endl;

/*

volScalarField infFilter(DiracInf*epsilon);

saaConcInfIntrpl *= infFilter;
saaConcInfIntrpl.correctBoundaryConditions();
DiracIntrpl *= infFilter;
DiracIntrpl.correctBoundaryConditions();
*/

for (int i=0; i<infIntrplRoundCell; i++)
{
	Info << "   Interpolation cell to vertices to cell: "<< i+1 << endl;
	const volPointInterpolation& pInterp = volPointInterpolation::New(mesh);
	const pointScalarField saaConcInfp =  pInterp.interpolate(saaConcInfIntrpl);
	const pointScalarField magLambdap =  pInterp.interpolate(DiracIntrpl);
	forAll(alpha1,celli)
	{
		saaConcInfIntrpl[celli] = interpolatePointToCell(saaConcInfp,celli);
		DiracIntrpl[celli] = interpolatePointToCell(magLambdap,celli);
	}
	saaConcInfIntrpl.correctBoundaryConditions();
	DiracIntrpl.correctBoundaryConditions();
}

volScalarField infFilter(DiracInf*epsilon);

saaConcInfIntrpl *= infFilter;
saaConcInfIntrpl.correctBoundaryConditions();
DiracIntrpl *= infFilter;
DiracIntrpl.correctBoundaryConditions();

saaConcInfIntrpl.oldTime() = saaConcInfIntrpl;
DiracIntrpl.oldTime() = DiracIntrpl;


surfaceTensorField NPft(mixture.nHatfv()*mixture.nHatfv());
surfaceTensorField Ift
(
	IOobject
	(
		"Ift",
		runTime.timeName(),
        mesh
	),
	mesh,
	tensor(1,0,0,0,1,0,0,0,1), //(XX, XY, XZ, YX, YY, YZ, ZX, ZY, ZZ)
	NPft.boundaryField().types()
);

//Info << "   Fake diffusion scaler: " << pow(cellWidthMin,2)/dt) << endl;

dimensionedScalar coefInfInrpl0 = (pow(cellWidthMin,2)/dt) * coefInfInrpl;
surfaceTensorField coefInfInrpl1((NPft.T()* coefInfInrpl0) & NPft);
surfaceScalarField coefInfInrpl2((fvc::interpolate(fvc::div(NPft.T() & mesh.Sf())) & (coefInfInrpl0*NPft)) & mesh.Sf());
volScalarField coefInfInrpl3(fvc::div(coefInfInrpl2));

for (int i=0; i<infIntrplRoundNormal; i++)
{
	Info << "   Interpolation in the normal direction: "<< i+1 << endl;

	fvScalarMatrix saaConcInfIntrplEqn
	(
		fvm::ddt(saaConcInfIntrpl)
		==
		fvm::laplacian(coefInfInrpl1,saaConcInfIntrpl)-fvm::div(coefInfInrpl2,saaConcInfIntrpl)//+fvm::Sp(coefInfInrpl3,saaConcInfIntrpl)
	);
	saaConcInfIntrplEqn.solve();
	fvScalarMatrix DiracIntrplEqn
	(
		fvm::ddt(DiracIntrpl)
		==
		fvm::laplacian(coefInfInrpl1,DiracIntrpl)-fvm::div(coefInfInrpl2,DiracIntrpl)//+fvm::Sp(coefInfInrpl3,DiracIntrpl)
	);
	DiracIntrplEqn.solve();
	
	saaConcInfIntrpl.oldTime() = saaConcInfIntrpl;
	DiracIntrpl.oldTime() = DiracIntrpl;
}


//saaConcInfArea = saaConcInfIntrpl / (DiracIntrpl + SMALL*dimensionedScalar("",dimensionSet(0,-1,0,0,0,0,0),1.0));
saaConcInfArea = saaConcInfIntrpl / (DiracIntrpl + SMALL/epsilon);

	
Info << "   Interface restriction" << endl;	
forAll(saaConcInfArea,celli)
{
	if (mag(cellDistInf[celli]) > (maxCellDist-double(1)+SMALL))
	{
		saaConcInfArea[celli] = double(0);	
	}
}

Info << "   Correction: internalField" << endl; // Limit cells value to have a possitive surface tension	

dimensionedScalar cstR("cstR",dimensionSet(1,2,-2,-1,-1,0,0),scalar(8.314));
dimensionedScalar cstT("cstT",dimensionSet(0,0,0,1,0,0,0),scalar(298.15));
dimensionedScalar cstVal = saaN * cstR * cstT * saaConcSurfMax;
//scalar SMALLmax = scalar(1e-2) * saaConcSurfMax.value();
//scalar cstLim = (1 - Foam::exp(-sigma0.value() / cstVal.value())) * saaConcSurfMax.value() - SMALLmax;
scalar sigmaMin = max(SMALL,sigmaMinLimit.value());
scalar cstLim = (1 - Foam::exp( (sigmaMin-sigma0.value()) / cstVal.value())) * saaConcSurfMax.value();
//Info << "	Correct cell values (lim=" << cstLim << ")" << endl;
//int nbrCorr=0;
forAll(saaConcInfArea,celli)
{
	if (saaConcInfArea[celli]>cstLim)
	{
		//++nbrCorr;
		saaConcInfArea[celli] = cstLim;	
	}
}
//Info << "		Nbr of corrections: " << nbrCorr << endl;

Info << "   Correction: boundaryField" << endl;
saaConcInfArea.correctBoundaryConditions();

}
