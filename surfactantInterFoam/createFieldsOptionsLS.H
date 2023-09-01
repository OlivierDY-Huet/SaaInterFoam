//Info<<"\n LS variables from solver dict" << endl;




Info<<"   Reading solverDict: " << "epsilon factor"<< endl;
const double epsilonFact = psiControls.lookupOrDefault<double>("epsilonFact", 1.4);
dimensionedScalar epsilon(epsilonFact*cellWidthMin);

Info<<"   Reading solverDict: " << "deltaTau factor"<< endl;
const double deltaTauFact = psiControls.lookupOrDefault<double>("deltaTauFact", 0.1);
dimensionedScalar deltaTau(deltaTauFact*cellWidthMin);

Info<<"   Reading solverDict: " << "relGradAlpha1LimWall"<< endl;
const double relGradAlpha1LimWall = psiControls.lookupOrDefault<double>("relGradAlpha1LimWall", 0.1);

Info<<"   Reading solverDict: " << "compLS"<< endl;
const double compLS = psiControls.lookupOrDefault<double>("compLS", 10);

Info<<"   Reading solverDict: " << "shCoef"<< endl;
const double shCoef = psiControls.lookupOrDefault<double>("shCoef", 0.5);

Info<<"\n Options for LS model" << endl;


const bool sorptionSwitch = transportProperties.lookupOrDefault("sorptionSwitch",true);
const bool saaInfIniSwitch = transportProperties.lookupOrDefault("saaInfIniSwitch",true);

/*
const int VOFLSSwitch = transportProperties.lookupOrDefault<int>("VOFLSSwitch", 3);

const int interpolateSSFField = transportProperties.lookupOrDefault<int>("interpolateSSFField", 0);

const int ReinitSwitch = transportProperties.lookupOrDefault<int>("ReinitSwitch", 1);

const bool advectPsiSwitch = transportProperties.lookupOrDefault("advectPsi",false); 
*/

Info<<"\n Options for Saa model" << endl;

const double saaWallEffCoef = transportProperties.lookupOrDefault<double>("saaWallEff", 1);
