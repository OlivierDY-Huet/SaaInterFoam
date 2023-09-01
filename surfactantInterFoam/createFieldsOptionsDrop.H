Info<<"\n Read initial condition for the drop" << endl;

dimensionedScalar dropRadius
(
	"dropRadius",
	dimensionSet(0,1,0,0,0,0,0),
	transportProperties
);

dimensionedVector dropCentre
(
	"dropCentre",
	dimensionSet(0,1,0,0,0,0,0),
	transportProperties
);

dimensionedVector dropVelocity
(
	"dropVelocity",
	dimensionSet(0,1,-1,0,0,0,0),
	transportProperties
);
