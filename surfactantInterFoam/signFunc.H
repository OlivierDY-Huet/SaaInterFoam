{
//https://www.eccomas2016.org/proceedings/pdf/9241.pdf

//signFunc = psi0/(pow(psi0*psi0 + pow(mag(gradPsi)+mag(gradPsiWall),2)*pow(cellWidth,2),0.5));
//signFunc = psi0/(pow(psi0*psi0 + pow(Gtot,2)*pow(cellWidth,2),0.5));	

volScalarField coefSign(max(mag(gradPsi)+mag(gradPsiWall),1.0));
volScalarField epsilonMod(coefSign*epsilon);
//volScalarField epsilonMod(max(Gtot,1.0)*epsilon);

signFunc = 	(psi/epsilonMod)+sin(M_PI*psi/epsilonMod)/M_PI;
forAll(psi,celli)
{
	if(psi[celli]<-epsilonMod[celli])
	{
		signFunc[celli] = -1.0;
	}
	else if(psi[celli]>epsilonMod[celli])
	{
		signFunc[celli] = 1.0;
	}
}
signFunc.correctBoundaryConditions();	
	
	
}
