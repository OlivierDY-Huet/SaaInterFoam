Info<< "\nExtract fields from interfaceProperties" << endl;
{
	
Info << "   Surface tension" << endl;
sigma = mixture.sigmaIP();
gradSigma = fvc::grad(sigma)*mixture.nearInterface1Epsilon();

//Info << "   nHat" << endl;
nHatInf = mixture.nHat();//fvc::reconstruct(mixture.nHatNWCf());

Info << "   Surface Tension force" << endl;
SurfaceTensionForce = mixture.surfaceTensionForce();

Info << "   Marangoni force" << endl;
MarangoniForce = coefMarangoni * mixture.MarangoniForce();

Info << "   Scaled Marangoni force" << endl;
MarangoniScaledForce = coefMarangoni * mixture.MarangoniDensityScaledForce();

Info << "   Scaled coefficient" << endl;
//densityScaleCoef = mixture.densityScaleCoef();
scaleCoefInf = mixture.Heaviside();

Info << "   Contact angle" << endl;

thetaDeg = mixture.thetaRad()*radToDeg();

Info << "   Curvature" << endl;

K = mixture.KCoef()*mixture.nearInterface1Epsilon();


gradHeaviside = fvc::grad(HeavisideInf);



}
