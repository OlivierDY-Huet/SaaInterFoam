Info<< "\nSolve U equation" << endl;


Info <<    "Correct boundary" << endl;
MRF.correctBoundaryVelocity(U);

Info << "   Relax" << endl;
fvVectorMatrix UEqn
(
	fvm::ddt(rho, U) 
  + fvm::div(rhoPhi, U)
  + MRF.DDt(rho, U)
  + turbulence->divDevRhoReff(rho, U)
 ==
	fvOptions(rho, U)
);
UEqn.relax();

Info << "   Constrain" << endl;
fvOptions.constrain(UEqn);


Info << "   Solve" << endl;
if (pimple.momentumPredictor())
{
	solve
	(
		UEqn
	 ==
		fvc::reconstruct
		(
			(
				mixture.surfaceTensionDensityScaledForcef() //mixture.surfaceTensionDensityScaledForcef()
			  + coefMarangoni * mixture.MarangoniDensityScaledForcef()
			  - ghf*fvc::snGrad(rho)
			  - fvc::snGrad(p_rgh)
			) * mesh.magSf()
		)
	);
	fvOptions.correct(U);
}
