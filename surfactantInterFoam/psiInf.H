Info<<"		psi seeding" << endl;
{

forAll(psi,celli)
{
	if (mag(psi[celli]/cellWidthMin.value()) < 0.5)//if (mag(psi[celli]) < epsilon.value()/2.0)
	{
		psi[celli] = psi0[celli];
	}
}

psi.correctBoundaryConditions();
}
