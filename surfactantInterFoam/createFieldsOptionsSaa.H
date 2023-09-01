//Info<<"\n LS variables from solver dict" << endl;


Info<<"   Reading transportProperties: " << "sigma0" << endl;
dimensionedScalar sigma0
("sigma0",dimensionSet(1,0,-2,0,0,0,0),readScalar(transportProperties.subDict("sigma").lookup("sigma0")));

Info<<"   Reading transportProperties: " << "sigmaMinLimit" << endl;
dimensionedScalar sigmaMinLimit("sigmaMinLimit",dimensionSet(1,0,-2,0,0,0,0),readScalar(transportProperties.subDict("sigma").lookup("sigmaMinLimit")));

Info<<"   Reading transportProperties: " << "saaN" << endl;
dimensionedScalar saaN("saaN",dimensionSet(0,0,0,0,0,0,0),readScalar(transportProperties.subDict("sigma").lookup("saaN")));



Info<<"   Reading transportProperties: " << "saaCMC" << endl;
dimensionedScalar saaCMC("saaCMC",dimensionSet(0,-3,0,0,1,0,0),transportProperties);

Info<<"   Reading transportProperties: " << "saaConcBulkIni" << endl;
dimensionedScalar saaConcBulkIni("saaConcBulkIni",dimensionSet(0,-3,0,0,1,0,0),transportProperties);

Info<<"   Reading transportProperties: " << "diffSaaBulk" << endl;
dimensionedScalar diffSaaBulk("diffSaaBulk",dimViscosity,transportProperties);

Info<<"   Reading transportProperties: " << "coefAd" << endl;
dimensionedScalar coefAd("coefAd",dimVelocity,transportProperties);

Info<<"   Reading transportProperties: " << "isothermLangmuir" << endl;
dimensionedScalar isothermLangmuir("isothermLangmuir",dimensionSet(0,3,0,0,-1,0,0),transportProperties);

Info<<"   Reading transportProperties: " << "saaConcSurfMax" << endl;
dimensionedScalar saaConcSurfMax("saaConcSurfMax",dimensionSet(0,-2,0,0,1,0,0),transportProperties);

Info<<"   Reading solverDict: " << "coefBarrier" << endl;
const double coefBarrier = transportProperties.lookupOrDefault<double>("coefBarrier", 0);

Info<<"   Reading solverDict: " << "coefDesorpBoost" << endl;
const double coefDesorpBoost = transportProperties.lookupOrDefault<double>("coefDesorpBoost", 1);




Info<<"   Reading solverDict: " << "coefBulkDiffSMALL" << endl;
const double coefBulkDiffSMALL = saaBulkControls.lookupOrDefault<double>("coefBulkDiffSMALL", 1e-3);

Info<<"   Reading solverDict: " << "coefBulkComp" << endl;
const double coefBulkComp = saaBulkControls.lookupOrDefault<double>("coefBulkComp", 1);

Info<<"   Reading solverDict: " << "coefBulkCorr" << endl;
const double coefBulkCorr = saaBulkControls.lookupOrDefault<double>("coefBulkCorr", 1);

Info<<"   Reading solverDict: " << "coefBulkCorrSMALL" << endl;
const double coefBulkCorrSMALL = saaBulkControls.lookupOrDefault<double>("coefBulkCorrSMALL", 1e-3);

Info<<"   Reading solverDict: " << "coefBulkCorrLim" << endl;
const double coefBulkCorrLim = saaBulkControls.lookupOrDefault<double>("coefBulkCorrLim", 1e2);



Info<<"   Reading solverDict: " << "coefInfComp" << endl;
const double coefInfComp = saaInfControls.lookupOrDefault<double>("coefInfComp", 1);

Info<<"   Reading solverDict: " << "coefInfCorr" << endl;
const double coefInfCorr = saaInfControls.lookupOrDefault<double>("coefInfCorr", 1);

Info<<"   Reading solverDict: " << "coefInfCorrSMALL" << endl;
const double coefInfCorrSMALL = saaInfControls.lookupOrDefault<double>("coefInfCorrSMALL", 1e-3);

Info<<"   Reading solverDict: " << "coefInfCorrLim" << endl;
const double coefInfCorrLim = saaInfControls.lookupOrDefault<double>("coefInfCorrLim", 1e2);



Info<<"   Reading solverDict: " << "coefMarangoni" << endl;
const double coefMarangoni = sigmaControls.lookupOrDefault<double>("coefMarangoni", 1);

Info<<"   Reading solverDict: " << "coefWallST" << endl;
const double coefWallST = sigmaControls.lookupOrDefault<double>("coefWallST", 0);

Info<<"   Reading solverDict: " << "infIntrplRoundCell" << endl;
const int infIntrplRoundCell = sigmaControls.lookupOrDefault<int>("infIntrplRoundCell", 1);

Info<<"   Reading solverDict: " << "infIntrplRoundNormal" << endl;
const int infIntrplRoundNormal = sigmaControls.lookupOrDefault<int>("infIntrplRoundNormal", 1);

Info<<"   Reading solverDict: " << "coefInfInrpl" << endl;
const double coefInfInrpl = sigmaControls.lookupOrDefault<double>("coefInfInrpl", 1);



