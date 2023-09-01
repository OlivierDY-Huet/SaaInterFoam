Info <<"Sink and source calculation"<<endl;

//Initialisation 
saaSourceSu = double(0)*saaSourceSu;
saaSourceSp = double(0)*saaSourceSp;
saaSinkSu = double(0)*saaSinkSu;
saaSinkSp = double(0)*saaSinkSp;

saaSourceSp1 = double(0)*saaSourceSp1;
saaSourceSp2 = double(0)*saaSourceSp2;

saaSourceSu.correctBoundaryConditions();
saaSourceSp.correctBoundaryConditions();
saaSinkSu.correctBoundaryConditions();
saaSinkSp.correctBoundaryConditions();
saaSourceSp1.correctBoundaryConditions();
saaSourceSp2.correctBoundaryConditions();


if (sorptionSwitch==true)
{
	
	//	//Processors
	Info <<"   Update processor values"<<endl;
	for (int iS=0; iS<sublayerIDs.size(); ++iS)
	{
		procSaaCBNew[thisProcNo][iS] = saaCBNew[sublayerIDs[iS]];
		/*
		if (procSaaCBNew[thisProcNo][iS]>saaCMC.value())
		{
			procSaaCBNew[thisProcNo][iS] = saaCMC.value();
		}
		*/
	}
	Pstream::gatherList(procSaaCBNew);
	Pstream::scatterList(procSaaCBNew);
	
	
	
	//	//Barrier
	
	/*
	volScalarField barrierBalance(saaConcInfArea/saaConcSurfMax);
	for (int i0=0; i0<infIDs.size(); ++i0) 
	{
		if (saaConcInfArea[infIDs[i0]]>saaConcSurfEq.value())
		{
			barrierBalance[infIDs[i0]] = max(saaConcSurfEq.value()/saaConcSurfMax.value()-coefDesorpBoost*(saaConcInfArea[infIDs[i0]]-saaConcSurfEq.value())/saaConcSurfMax.value(),double(0.0));
		}
	}
	barrierBalance.correctBoundaryConditions();
	sorptBarrierInf = mixture.nearInterface1Epsilon() * Foam::exp(-coefBarrier*barrierBalance);
	*/
	
	volScalarField barrierBalance(saaConcInfArea/saaConcSurfMax);
	for (int i0=0; i0<infIDs.size(); ++i0) 
	{
		if (saaConcInfArea[infIDs[i0]]>saaConcSurfEqMax.value())
		{
			barrierBalance[infIDs[i0]] = saaConcSurfEqMax.value()/saaConcSurfMax.value()-coefDesorpBoost*(saaConcInfArea[infIDs[i0]]-saaConcSurfEqMax.value())/saaConcSurfMax.value();
		}
	}
	barrierBalance.correctBoundaryConditions();
	sorptBarrierInf = mixture.nearInterface1Epsilon() * Foam::exp(-coefBarrier*max(barrierBalance,double(0.0)));
	
	//Source and sink calculation
	Info <<"   Sink/Source exchange"<<endl;
	{
	int neighSublayerIdx = -1;
	for (int i0=0; i0<infIDs.size(); ++i0) 
	{

		scalar coefAdVal = coefAd.value();
		scalar coefDeVal = coefAd.value()/(isothermLangmuir.value()*saaConcSurfMax.value());
		
		//Info <<"	Interface cell ID:"<< infIDs[i0] <<endl;
		for (int kS=0; kS<infLinkIdx[i0].size(); ++kS)
		{
			
			
			if (infLinkIdx[i0][kS] != -1) //-1: no pairing
			{
				
				int SPOri = infLinkOri[i0][kS]; // Sublayer Processor Origin
				int SPIdx = infLinkIdx[i0][kS]; // Sublayer Processor Index
				
				scalar w = infLinkWeight[i0][kS]; //weight of the operation
				scalar rVol = meshV[infIDs[i0]]/procSublayerVol[SPOri][SPIdx]; //Ratio of volume
				
				scalar CBLim = min(procSaaCBNew[SPOri][SPIdx],saaCMC.value());
				scalar CBExc = min(procSaaCBNew[SPOri][SPIdx]-saaCMC.value(),0);
						
				saaSourceSu[infIDs[i0]] += w*sorptBarrierInf[infIDs[i0]]*(coefAdVal*DiracInf[infIDs[i0]]*CBLim);
				saaSourceSp[infIDs[i0]] -= w*sorptBarrierInf[infIDs[i0]]*(coefAdVal*CBLim/saaConcSurfMax.value()+coefDeVal);
				
				saaSourceSp1[infIDs[i0]] -= w*sorptBarrierInf[infIDs[i0]]*(coefAdVal*CBLim/saaConcSurfMax.value());
				saaSourceSp2[infIDs[i0]] -= w*sorptBarrierInf[infIDs[i0]]*coefDeVal;
				
				//scalar sinkSuVal = rVol*w*sorptBarrierInf[infIDs[i0]]*(coefDeVal*saaCINew[infIDs[i0]]);
				scalar sinkSuVal = rVol*w*sorptBarrierInf[infIDs[i0]]*(coefDeVal*saaCINew[infIDs[i0]]+(coefAdVal*(DiracInf[infIDs[i0]]-saaCINew[infIDs[i0]]/saaConcSurfMax.value()))*CBExc);
				scalar sinkSpVal = -rVol*w*sorptBarrierInf[infIDs[i0]]*(coefAdVal*(DiracInf[infIDs[i0]]-saaCINew[infIDs[i0]]/saaConcSurfMax.value()));
				
				if (SPOri == thisProcNo)
				{
					saaSinkSu[sublayerIDs[SPIdx]] += sinkSuVal;
					saaSinkSp[sublayerIDs[SPIdx]] += sinkSpVal;
				}
				else
				{
					neighSublayerIdx += 1;
					procSaaSinkSuCorr[thisProcNo][neighSublayerIdx] = sinkSuVal;
					procSaaSinkSpCorr[thisProcNo][neighSublayerIdx] = sinkSpVal;
				}				
			}
		}
	}
	}
	
	saaSourceSu.correctBoundaryConditions();
	saaSourceSp.correctBoundaryConditions();
	saaSinkSu.correctBoundaryConditions();
	saaSinkSp.correctBoundaryConditions();
	saaSourceSp1.correctBoundaryConditions();
	saaSourceSp2.correctBoundaryConditions();
	
	
	//Sink correction
	Info <<"   Sink correction" <<endl;
	Pstream::gatherList(procSaaSinkSuCorr);
	Pstream::gatherList(procSaaSinkSpCorr);
	Pstream::scatterList(procSaaSinkSuCorr);
	Pstream::scatterList(procSaaSinkSpCorr);
	for(int pP=0; pP<procPatchIDs.size(); ++pP)
	{
		int refp=procPatchRefNo[thisProcNo][pP];
		if (procSaaSinkSuCorr[refp].size()>0)
		{
			for (int nCorr=0; nCorr<procSaaSinkSuCorr[refp].size(); ++nCorr)
			{
				if (procSublayerProcOriCorr[refp][nCorr] == thisProcNo)
				{
					int idxCorr = procSublayerIdxCorr[refp][nCorr];
					saaSinkSu[sublayerIDs[idxCorr]] += procSaaSinkSuCorr[refp][nCorr];
					saaSinkSp[sublayerIDs[idxCorr]] += procSaaSinkSpCorr[refp][nCorr];
				}
			}
		}
	}
	saaSinkSu.correctBoundaryConditions();
	saaSinkSp.correctBoundaryConditions();

}


Info <<"Triple line"<<endl;

saaWallSp = double(0)*saaWallSp; 
if (saaWallEffCoef>SMALL)
{	
	const fvPatchList& patches = mesh.boundary();
    forAll(patches, patchi)
    {
        const fvPatch& curPatch = patches[patchi];
        if (isType<wallFvPatch>(curPatch))
        {	
			Info <<"	Wall patch: initialisation"<<endl;
			const fvPatchVectorField& Up = U.boundaryField()[patchi];
			vectorField nf = curPatch.nf();
				
			Info <<"	Wall patch: effective velocity"<<endl;
			vectorField Uwall = Up.patchInternalField() - Up;
			scalarField uWallEff = (nf & Uwall);
			
			Info <<"	Wall patch: saa "<<endl;
			forAll(curPatch,facei)
			{
				if (uWallEff[facei]>0)
				{
					label faceCelli = curPatch.faceCells()[facei];
					//scalar DiracWall = double(1)/(double(2)*cellWidth[faceCelli]);
						
					saaWallSp[faceCelli] += saaWallEffCoef*uWallEff[facei]*DiracWall[faceCelli];
					
				}
			}			
        }
    }
}
