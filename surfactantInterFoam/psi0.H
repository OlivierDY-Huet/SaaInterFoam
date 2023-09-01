Info<<"		psi0" << endl;
{

gradAlpha1= fvc::grad(alpha1);

//Raeini, A.Q.; Blunt, M.J.; Bijeljic, B. Modelling two-phase flow in porous media at the pore scale using the volume-of-fluid method. J. Comput. Phys. 2012, 231, 5653–5668.

/*
alpha1Sharp = (1.0/(1.0-shCoef))*(min(max(alpha1,shCoef/2.0),1.0-shCoef/2.0)-shCoef/2.0);
alpha1Sharp.correctBoundaryConditions();
*/

psi0 = (double(2.0)*alpha1-double(1.0))*(epsilon/2);
psi0.correctBoundaryConditions();




}
