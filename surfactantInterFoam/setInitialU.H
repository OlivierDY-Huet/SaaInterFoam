{
Info<<"		Reset U" << endl;

#include "cellWidth.H"
#include "meshLS.H"

/*
dimensionedScalar centreX = dropCentre[0];
dimensionedScalar centreY = dropCentre[1];
dimensionedScalar centreZ = dropCentre[2];


scalarField dX = X_Coor.field() - centreX.value();
scalarField dY = Y_Coor.field() - centreY.value();
scalarField dZ = Z_Coor.field() - centreZ.value();
scalarField dR = (dropRadius.value() + epsilon.value()) - pow(pow(dX,2) + pow(dY,2) + pow(dZ,2),0.5);
*/
volVectorField dXYZ(XYZ_Coor-dropCentre);
volScalarField dR(dropRadius-mag(dXYZ));

forAll(U,celli)
{
	if(dR[celli] < -epsilon.value())
		U[celli] = double(0.0)*dropVelocity.value();
	else if(dR[celli] > epsilon.value())
		U[celli] = double(1.0)*dropVelocity.value();
	else
		U[celli] = (double(1.0)/double(2.0))*(double(1.0)+(dR[celli]/epsilon.value())+Foam::sin(M_PI*dR[celli]/epsilon.value())/M_PI)*dropVelocity.value();
}

U.correctBoundaryConditions();

}
