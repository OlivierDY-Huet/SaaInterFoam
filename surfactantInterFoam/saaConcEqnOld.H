Info<< "\nSolve surfactant equations" << endl;
{


surfaceScalarField coefBulkH(fvc::interpolate(1.0-HeavisideInf));
surfaceScalarField coefInfH(fvc::interpolate((1.0-epsilon*DiracInf)*sign(-signFunc))); //surfaceScalarField coefInfH(fvc::interpolate(1.0-2.0*HeavisideInf));

//Schemes

//surfaceScalarField phir(phic_*(mixture.nHatfv() & mesh.Sf()));
phirField=fvc::reconstruct(phir_); 

surfaceScalarField phirBulk(coefBulkComp*phir_*coefBulkH);
surfaceScalarField phirInf(coefInfComp*phir_*coefInfH);

cmpCB = fvc::reconstruct(phirBulk);
cmpCI = fvc::reconstruct(phirInf);

//advection-diffusion
Info << "   Advection-diffusion" << endl;

surfaceScalarField phisBulk = (diffSaaBulk*fvc::snGrad(HeavisideInf)/fvc::interpolate(mag(HeavisideInf) + coefBulkDiffSMALL))*mesh.magSf();
phisBulkField = fvc::reconstruct(phisBulk);
//zeroField divU;

Info << "   	Advection-diffusion: c_B" << endl;
fvScalarMatrix saaConcBulkADEqn
(
	fvm::ddt(saaConcBulk)
	+ fvm::div(phi,saaConcBulk)
	+ fvm::div(phirBulk,saaConcBulk)	
	==
	fvm::laplacian(diffSaaBulk,saaConcBulk)
	- fvm::div(phisBulk,saaConcBulk)
);
saaConcBulkADEqn.solve();

Info << "   	Advection-diffusion: c_I" << endl;
fvScalarMatrix saaConcInfADEqn
(
	fvm::ddt(saaConcInf)
	+ fvm::div(phi,saaConcInf)
	+ fvm::div(phirInf,saaConcInf)
);
saaConcInfADEqn.solve();


saaCBNew = saaConcBulk; 
saaCINew = saaConcInf;
saaConcBulk = saaConcBulk.oldTime();
saaConcInf = saaConcInf.oldTime();
saaConcBulk.oldTime()=saaCBNew;
saaConcInf.oldTime()=saaCINew;
saaConcBulk = saaCBNew;
saaConcInf = saaCINew;

Info << "   Correction" << endl;


surfaceScalarField innerInff(fvc::interpolate(innerInf));
innerInff = pos0(innerInff-double(0.9));

//volTensorField NPt(mixture.nHat()*mixture.nHat());
surfaceTensorField NPft(mixture.nHatHeavisidefv()*mixture.nHatHeavisidefv());//surfaceTensorField NPft(mixture.nHatfv()*mixture.nHatfv()); 

surfaceScalarField coefBulkGradf(cellWidthMin * innerInff *fvc::interpolate(1.0/(HeavisideInf+coefBulkCorrSMALL)) *fvc::snGrad(HeavisideInf) *  mesh.magSf());
surfaceScalarField coefInfGradf(cellWidthMin * innerInff *fvc::interpolate(1.0/(DiracInf+coefInfCorrSMALL/epsilon)) *fvc::snGrad(DiracInf) *  mesh.magSf());

coefBulkGrad = fvc::reconstruct(coefBulkGradf);
coefInfGrad =  fvc::reconstruct(coefInfGradf);

surfaceScalarField coefBulkLimf(fvc::interpolate((-signFunc)*(double(1)-pos0(HeavisideInf-SMALL))) * mixture.nHatfv() & mesh.Sf());
surfaceScalarField coefInfLimf(fvc::interpolate((-signFunc)*(double(1)-pos0(DiracInf*epsilon-SMALL))) * mixture.nHatfv() & mesh.Sf());

coefBulkLim = fvc::reconstruct(coefBulkLimf);
coefInfLim =  fvc::reconstruct(coefInfLimf);

int corrNbr = 1;//int(10);
double corrNrrInv = 1.0/double(corrNbr);


surfaceTensorField bulkCorrLapf(corrNrrInv * innerInff *  coefBulkCorr * phic_ * cellWidthMin * NPft);
surfaceScalarField bulkCorrGradf(corrNrrInv * phic_  * (coefBulkCorr * coefBulkGradf + coefBulkCorrLim* coefBulkLimf));

surfaceTensorField infCorrLapf(corrNrrInv * innerInff * coefInfCorr * phic_ * cellWidthMin * NPft);
surfaceScalarField infCorrGradf(corrNrrInv * phic_  * (coefInfCorr * coefInfGradf + coefInfCorrLim * coefInfLimf));


for (int corrLoopCount = 0; corrLoopCount < corrNbr; corrLoopCount++) 
{
	Info << "   Correction: c_B" << endl;
	fvScalarMatrix saaConcBulkCorrEqn
	(
		fvm::ddt(saaConcBulk)
		==
		fvm::laplacian(bulkCorrLapf,saaConcBulk) - fvm::div(bulkCorrGradf,saaConcBulk)
	);
	saaConcBulkCorrEqn.solve();

	Info << "   Correction: c_I" << endl;
	fvScalarMatrix saaConcInfCorrEqn
	(
		fvm::ddt(saaConcInf)
		==
		fvm::laplacian(infCorrLapf,saaConcInf) - fvm::div(infCorrGradf,saaConcInf)
	);
	saaConcInfCorrEqn.solve();

	saaCBNew = saaConcBulk; 
	saaCINew = saaConcInf;
	saaConcBulk = saaConcBulk.oldTime();
	saaConcInf = saaConcInf.oldTime();
	saaConcBulk.oldTime()=saaCBNew;
	saaConcInf.oldTime()=saaCINew;
	saaConcBulk = saaCBNew;
	saaConcInf = saaCINew;
}


//sink-source + corrections
Info << "   Sink/source loop" << endl;

#include "cellInfDist.H"
#include "saaSourceSinkLink.H"
#include "saaSourceSink.H"

scalar saaCBErrMax = double(1000);
scalar saaCIErrMax = double(1000);
scalar maxErr = double(1000);
label saaConvLoopMax(runTime.controlDict().get<label>("saaConvLoopMax"));
scalar saaConvLoopTol(runTime.controlDict().get<scalar>("saaConvLoopTol"));

for (int saaLoopCount = 0; saaLoopCount < saaConvLoopMax; saaLoopCount++) 
{
	Info << "   Sink/source loop: c_B" << endl;
	fvScalarMatrix saaConcBulkSSEqn
	(
		fvm::ddt(saaConcBulk)
		==
		fvm::Sp(saaSinkSp,saaConcBulk) + saaSinkSu
	);
	saaConcBulkSSEqn.solve();
	
	Info << "   Sink/source loop: c_I" << endl;
	fvScalarMatrix saaConcInfSSEqn
	(
		fvm::ddt(saaConcInf)
		==
		saaSourceSu + fvm::Sp(saaSourceSp,saaConcInf)
		- fvm::Sp(saaWallSp,saaConcInf)
	);
	saaConcInfSSEqn.solve();
	
	Info << "   Sink/source loop: c_W" << endl;
	fvScalarMatrix saaConcWallSSEqn
	(
		fvm::ddt(saaConcWall)
		==
		saaWallSp * saaCINew
	);
	saaConcWallSSEqn.solve();
	
	
	saaCBErrMax = max(mag(saaConcBulk.field() - saaCBNew.field())/(mag(saaCBNew.field())+SMALL));
	saaCIErrMax = max(mag(saaConcInf.field() - saaCINew.field())/(mag(saaCINew.field())+SMALL));
	reduce(saaCBErrMax, maxOp<scalar>());
	reduce(saaCIErrMax, maxOp<scalar>());
	maxErr = max(saaCBErrMax,saaCIErrMax);
	Info << "   Sink/source loop: Maximum error =" << maxErr << endl;

	Info << "   Sink/source loop: Update concentrations for source/sink" << endl;
	saaCBNew = saaConcBulk; 
	saaCINew = saaConcInf;
	saaCWNew = saaConcWall;
	#include "saaSourceSink.H"
	
	#include "saaConcInfArea.H" ///////////////
	
	Info << "   Sink/source loop: reset saa Fields" << endl; 
	saaConcBulk = saaConcBulk.oldTime();
	saaConcInf = saaConcInf.oldTime();
	saaConcWall = saaConcWall.oldTime();
	
	if (maxErr<saaConvLoopTol)
	{
		break;
	}
}

saaConcBulk = saaCBNew; 
saaConcInf = saaCINew;
saaConcWall = saaCWNew;

Info << "   Convergence loop: final Source/Sink" << endl;
saaSourceSpC1 = saaSourceSp1 * saaConcInf;
saaSourceSpC2 = saaSourceSp2 * saaConcInf;
saaSinkSpC = saaSinkSp * saaConcBulk;
saaSource = saaSourceSu + saaSourceSp * saaConcInf; 
saaSink = saaSinkSp * saaConcBulk + saaSinkSu;
saaWallSpC = saaWallSp * saaConcInf;




}
