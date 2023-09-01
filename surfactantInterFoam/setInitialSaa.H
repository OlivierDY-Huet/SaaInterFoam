Info<<"\nInitialise saa for bulk and interface" << endl;

{

Info<<"   Initialise saa bulk" << endl;

#include "Dirac.H"

saaConcBulk = saaConcBulkIni * HeavisideInf;//saaConcBulk = saaConcBulkIni * alpha1;
//saaConcBulk = saaConcBulkIni +  saaConcBulkIni * (HeavisideInf-1)/(HeavisideInf+SMALL);
forAll(saaConcBulk,celli)
{
	if(saaConcBulk[celli]<double(0.0))
	{
		saaConcBulk[celli] = double(0.0);
	}
}
saaConcBulk.correctBoundaryConditions();

Info<<"   Initialise saa interface" << endl;
//Equilibrium
if (coefAd.value() > 0.0)
{
	saaConcSurfEq  = saaConcBulkIni*saaConcSurfMax/(saaConcBulkIni+(double(1)/isothermLangmuir));	
}
else
{
	saaConcSurfEq= (1/2)*saaConcSurfMax; //Defalt value that wont't be used anyway
}
saaConcSurfEqMax  = saaCMC*saaConcSurfMax/(saaCMC+(double(1)/isothermLangmuir));

//Choose value
if (saaInfIniSwitch == false)
{
	Info << "   Initialise saa interface: Not Equilibrium value" << endl;

}
else
{
	Info << "   Initialise saa interface: Equilibrium value" << endl;
	saaConcSurfIni = saaConcSurfEq;
}
Info << "   Initialise saa interface: " << saaConcSurfIni.value() << " * DiracInf "<< endl;

//Set value
saaConcInf = saaConcSurfIni * DiracInf;
saaConcInf.correctBoundaryConditions();

saaCBNew = saaConcBulk; 
saaCINew = saaConcInf;


}
