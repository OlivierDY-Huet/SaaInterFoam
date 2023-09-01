Info<<"\nCorrect LS" << endl;


//Info<<"   Initialisation of variables" << endl;
#include "psi0.H"

scalar loopNbr = epsilonFact/deltaTauFact;
Info<<"   Initialisation of psi with psi0" << LSIniSwitch <<endl;
if(LSIniSwitch==true)
{
	psi = psi0;
	psi.correctBoundaryConditions();
}
else
{
	#include "psiSeeding.H"
	loopNbr = (loopNbr/2.0);
}


Info<<"   Loop" << endl;
volScalarField psiOld = psi;


scalar psiMax = double(0.0);
scalar psiRelErr = double(1000);
scalar gradPsiErr = double(1000);

int count = 0;
for (int LSLoopCount=0; LSLoopCount<int(loopNbr); LSLoopCount++) // for (int LSLoopCount = 0; LSLoopCount < LSConvLoopMax; LSLoopCount++) //
{

	++count;
	
	//gradPsi
	gradPsi = fvc::grad(psi);	
	#include "gradPsiWall.H"
	
	//Gtot
	#include "Gtot.H"
	
	//Sign function
	#include "signFunc.H"
	

	//psi calculation
	//psi = psi - coefBC * deltaTau * signFunc * (pow(Gtot,2)-1.0);
	//psi = psi - coefBC * deltaTau * signFunc * (Gtot-1.0);
	//psi = psi - coefG * deltaTau * signFunc * (Gtot-1.0);
	//psi = psi - deltaTau * signFunc * (pow(Gtot,2)-1.0);
	//psi = psi - deltaTau * signFunc * (Gtot-1.0);
	psi = psi - deltaTau * signFunc * mag(signFunc) * (pow(mag(Gtotv),2)-1.0);
	psi.correctBoundaryConditions();

	//Error
	volScalarField checkErrorRelPsi = pos0((epsilon - mag(psi))/epsilon);
	volScalarField checkErrorGradPsi = pos0((0.5*epsilon - mag(psi))/(0.5*epsilon));

	psiRelErr = gMax(checkErrorRelPsi.field()*(mag(psi.field() - psiOld.field()))/(mag(psiOld.field())+VSMALL));
	//gradPsiErr = gMax(checkErrorGradPsi.field()*mag(mag(gradPsi)+mag(gradPsiWall)-1.0));
	gradPsiErr = gMax(checkErrorGradPsi.field()*mag(mag(Gtotv)-1.0));
	
	//Next round
	psiOld = psi;
}

Info<<nl<< "   Round " << count <<": max relative psi error = "<< psiRelErr << " ; max Gtot error "<< gradPsiErr << endl;
Info<<"\nCorrect LS" << endl;

psiInf = pos0(epsilon-mag(psi))*psi+(1-pos0(epsilon-mag(psi)))*sign(psi)*epsilon;
psiInf.correctBoundaryConditions();

psiRel = psi/cellWidthMin;
psiRel.correctBoundaryConditions();

LSIniSwitch = false;
