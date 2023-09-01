Info <<"Source-Sink link"<<endl;


Info <<"   References for processors"<<endl;

int totProcNo = Pstream::nProcs();
int thisProcNo = Pstream::myProcNo();


const fvPatchList& patches = mesh.boundary();
int nProcPatch=0;
int nProcElements=0;
std::vector<int> procPatchIDs;
forAll(patches, patchi)
{
	const fvPatch& curPatch = patches[patchi];
	if (isType<processorFvPatch>(curPatch))
	{
		nProcPatch += 1;
		nProcElements += curPatch.size();
		procPatchIDs.push_back(patchi); 
	}
}

List< List<int> > procPatchRefNo(totProcNo); //Processor references
List<int> procPatchRefNoP(nProcPatch); 
for(int pP=0; pP<procPatchIDs.size(); ++pP)
{
	const fvPatch& curPatch = patches[procPatchIDs[pP]];
	// Processor reference of the patch
	word patchName = curPatch.name();
	//Info << "Processor patch name:" << patchName << endl; //Template: ["procBoundary",thisProcNo,"to",patchProcessorNo]
	//const processorPolyPatch& procPolyPatch = refCast<const processorPolyPatch>(patches[procPatchIDs[pP]]); //https://www.cfd-online.com/Forums/openfoam-programming-development/136424-neighbouring-processor-number-patch.html
	//procPatchRefNoP[pP] = procPolyPatch.neighbProcNo();
	int strEnd = patchName.find_last_of("to");
	std::string patchProcessorNoStr = patchName.substr(strEnd+1);
	int patchProcessorNo = std::stoi(patchProcessorNoStr);
	procPatchRefNoP[pP] = patchProcessorNo;
}
procPatchRefNo[thisProcNo]=procPatchRefNoP;

//	//Gather - Scatter
Pstream::gatherList(procPatchRefNo);
Pstream::scatterList(procPatchRefNo);


Info <<"   Find interface"<<endl;

/*
std::vector<int> infIDs;
for (meshLvl=int(0); meshLvl<int(1); ++meshLvl)
{
	cellDistInfLvl=cellDistInfList[meshLvl];
	for (std::list<int>::iterator itS=cellDistInfLvl.begin(); itS!=cellDistInfLvl.end(); ++itS)
	{
			infIDs.push_back(*itS);
	}	
}
*/
std::vector<int> infIDs;
cellDistInfLvl=cellDistInfList[0];
for (std::list<int>::iterator itS=cellDistInfLvl.begin(); itS!=cellDistInfLvl.end(); ++itS)
{
	infIDs.push_back(*itS);
}	


Info <<"   Find sublayer "<<endl;
std::vector<int> sublayerIDs;
for (meshLvl=int(1); meshLvl<int(5); ++meshLvl)
{
	cellDistInfLvl=cellDistInfList[meshLvl];
	for (std::list<int>::iterator itS=cellDistInfLvl.begin(); itS!=cellDistInfLvl.end(); ++itS)
	{
		if (alpha1[*itS]>lim1lpha)
		{
			sublayerIDs.push_back(*itS);
		}
	}	
}


Info <<"   Initialisation "<<endl;

const scalarField& meshV = mesh.cellVolumes();

//Global
List< List<scalar> > procSublayerVol(totProcNo);
List< List<point> > procSublayerPosXYZ(totProcNo);
//Local
List<int> procSublayerIdxP(sublayerIDs.size());
List<scalar> procSublayerVolP(sublayerIDs.size());
List<point> procSublayerPosXYZP(sublayerIDs.size());
//Fuse: local + neighboors
std::vector<int> procSublayerProcOriF;
std::vector<int> procSublayerIdxF;

// Fill local and fused variables
for (int iS=0; iS<sublayerIDs.size(); ++iS)
{
	//Local
	procSublayerIdxP[iS] = iS;
	procSublayerVolP[iS]=meshV[sublayerIDs[iS]];
	procSublayerPosXYZP[iS]=mesh.C().internalField()[sublayerIDs[iS]];
	//Fuse
	procSublayerProcOriF.push_back(thisProcNo);
	procSublayerIdxF.push_back(iS);
}

// Transfer globally
procSublayerVol[thisProcNo] = procSublayerVolP;
procSublayerPosXYZ[thisProcNo] = procSublayerPosXYZP;
Pstream::gatherList(procSublayerVol);
Pstream::gatherList(procSublayerPosXYZ);
Pstream::scatterList(procSublayerVol);
Pstream::scatterList(procSublayerPosXYZ);

// Fuse the local sublayer to neighboor sublayers
for(int pP=0; pP<procPatchIDs.size(); ++pP)
{
	int refp=procPatchRefNo[thisProcNo][pP];
	for(int nF=0; nF<procSublayerVol[refp].size(); ++nF)
	{
		procSublayerProcOriF.push_back(refp);
		procSublayerIdxF.push_back(nF);
	}
}

Info <<"   Link "<<endl;
//Link between the interface and the sublayer  
std::vector<std::vector<int> > infLinkOri;
std::vector<int> list1Ori;
std::vector<int> list2Ori;
std::vector<int> listLinkOri;
std::vector<std::vector<int> > infLinkIdx;
std::vector<int> list1Idx;
std::vector<int> list2Idx;
std::vector<int> listLinkIdx;
std::vector<std::vector<scalar> > infLinkWeight;

std::vector<scalar> list1Dist;
std::vector<scalar> list2Dist; 
std::vector<scalar> listLinkDist;

std::vector<int> sublayerProcOriCorr;
std::vector<int> sublayerProcIdxCorr;

//scalar errLinkInf(runTime.controlDict().get<scalar>("errLinkInf"));
for (int i0=0; i0<infIDs.size(); ++i0)
{
	list1Ori.clear();
	list2Ori.clear();
	listLinkOri.clear();
	list1Idx.clear();
	list2Idx.clear();
	listLinkIdx.clear();
	list1Dist.clear();
	list2Dist.clear();

	scalar lim1 = cellWidth[infIDs[i0]]*double(10);//avgCellWidth.value()*5;
	scalar lim2 = lim1;
	
	scalar errLinkInf = cellWidth[infIDs[i0]]/double(10);
	point pt1 = mesh.C().internalField()[infIDs[i0]];
	for (int jS=0; jS<procSublayerIdxF.size(); ++jS)
	{
		point pt2 = procSublayerPosXYZ[procSublayerProcOriF[jS]][procSublayerIdxF[jS]];
		scalar distR = mag(pt2-pt1);

		if (distR<=lim1+errLinkInf)
		{
			if (distR<lim1-errLinkInf)
			{
				lim2=lim1;
				list2Ori=list1Ori;
				list1Ori.clear();
				list2Idx=list1Idx;
				list1Idx.clear();
				lim1=distR;
				list2Dist = list1Dist;
				list1Dist.clear();  
			}
			list1Ori.push_back(procSublayerProcOriF[jS]);
			list1Idx.push_back(procSublayerIdxF[jS]);
			list1Dist.push_back(distR);
		}
		else if (distR<=lim2+errLinkInf)
		{
			if (distR<lim2-errLinkInf)
			{
				list2Ori.clear();
				list2Idx.clear();
				lim2=distR;
				list2Dist.clear();
			}
			list2Ori.push_back(procSublayerProcOriF[jS]);
			list2Idx.push_back(procSublayerIdxF[jS]);	
			list2Dist.push_back(distR);
		}
	}
	//Fuse lists A and B
	listLinkOri = list1Ori;
	listLinkOri.insert(listLinkOri.end(),list2Ori.begin(),list2Ori.end());
	listLinkIdx = list1Idx;
	listLinkIdx.insert(listLinkIdx.end(),list2Idx.begin(),list2Idx.end());
	listLinkDist = list1Dist;
	listLinkDist.insert(listLinkDist.end(),list2Dist.begin(),list2Dist.end());
	if (listLinkOri.empty())
	{
		listLinkOri.push_back(label(-1));
		listLinkIdx.push_back(label(-1));
	}
	else
	{
		for (int kS=0; kS<listLinkIdx.size(); ++kS)
		{
			if(listLinkOri[kS] != thisProcNo)
			{
				sublayerProcOriCorr.push_back(listLinkOri[kS]);
				sublayerProcIdxCorr.push_back(listLinkIdx[kS]);
			}
		}
	}
	scalar w = 1/scalar(listLinkOri.size());
	//Info <<"	Size :" << listLinkOri.size() <<" ; Weight:" << w <<endl; /////////////////////////////////////
	std::vector<scalar> listLinkWeight(listLinkOri.size(),w);
	infLinkOri.push_back(listLinkOri);
	infLinkIdx.push_back(listLinkIdx);
	infLinkWeight.push_back(listLinkWeight);
}

Info <<"   Initialisation of sub terms"<<endl;
//Initialise processor lists 
//	//for correction
List< List<scalar> > procSublayerProcOriCorr(totProcNo);
List< List<scalar> > procSublayerIdxCorr(totProcNo);
List< List<scalar> > procSaaSinkSuCorr(totProcNo);
List< List<scalar> > procSaaSinkSpCorr(totProcNo);
List< List<scalar> > procSublayerLinkCorr(totProcNo);

List<scalar> procSublayerProcOriCorrP(sublayerProcIdxCorr.size());
List<scalar> procSublayerIdxCorrP(sublayerProcIdxCorr.size());
List<scalar> procSaaSinkSuCorrP(sublayerProcIdxCorr.size());
List<scalar> procSaaSinkSpCorrP(sublayerProcIdxCorr.size());
List<scalar> procSublayerLinkCorrP(sublayerProcIdxCorr.size());

for (int nCorr=0; nCorr<sublayerProcIdxCorr.size(); ++nCorr)
{
	procSublayerProcOriCorrP[nCorr] = sublayerProcOriCorr[nCorr];
	procSublayerIdxCorrP[nCorr] = sublayerProcIdxCorr[nCorr];
	//Info <<"	Transfert for corr: Ori: " << procSublayerProcOriCorrP[nCorr]  << " ; Idx" << procSublayerIdxCorrP[nCorr] <<endl; /////////////////////////////////////
}

procSublayerProcOriCorr[thisProcNo] = procSublayerProcOriCorrP;
procSublayerIdxCorr[thisProcNo] = procSublayerIdxCorrP;
procSaaSinkSuCorr[thisProcNo] = procSaaSinkSuCorrP;
procSaaSinkSpCorr[thisProcNo] = procSaaSinkSpCorrP;
procSublayerLinkCorr[thisProcNo] = procSublayerLinkCorrP;

//	//for source/sink
List< List<scalar> > procSaaCBNew(totProcNo);
List<scalar> procSaaCBNewP(sublayerIDs.size());
//List< List<scalar> > procSaaCoefAd(totProcNo);
List<scalar> procSaaCoefAdP(sublayerIDs.size());
procSaaCBNew[thisProcNo] = procSaaCBNewP;
//procSaaCoefAd[thisProcNo] = procSaaCoefAdP;

//	//Gather-scatter
Pstream::gatherList(procSublayerProcOriCorr);
Pstream::gatherList(procSublayerIdxCorr);
Pstream::gatherList(procSaaSinkSuCorr);
Pstream::gatherList(procSaaSinkSpCorr);

Pstream::gatherList(procSublayerLinkCorr);
Pstream::gatherList(procSaaCBNew);
//Pstream::gatherList(procSaaCoefAd);

Pstream::scatterList(procSublayerProcOriCorr);
Pstream::scatterList(procSublayerIdxCorr);
Pstream::scatterList(procSaaSinkSuCorr);
Pstream::scatterList(procSaaSinkSpCorr);
Pstream::scatterList(procSublayerLinkCorr);
Pstream::scatterList(procSaaCBNew);
//Pstream::scatterList(procSaaCoefAd);

Info <<"   Link evaluation "<<endl;
//	//To local + store neighbour values
checkInfLink = 0 * checkInfLink;
int neighSublayerIdx = -1;
//Info <<"Corr list size :" << procSublayerLinkCorr[thisProcNo].size() <<endl; /////////////////////////////////////
for (int i0=0; i0<infLinkIdx.size(); ++i0) 
{
	//Info <<"New inf cell" <<endl;
	for (int kS=0; kS<infLinkIdx[i0].size(); ++kS)
	{
		if (infLinkIdx[i0][kS] != -1) //-1: no pairing
		{
			int SPOri = infLinkOri[i0][kS]; // Sublayer Processor Origin
			int SPIdx = infLinkIdx[i0][kS]; // Sublayer Processor Index
			if (SPOri == thisProcNo)
			{
				checkInfLink[sublayerIDs[SPIdx]] += infLinkWeight[i0][kS];
			}
			else //else if(SPOri != -1)
			{
				neighSublayerIdx += 1;
				//Info <<"	n " << neighSublayerIdx <<" : Ori=" << SPOri <<" ; Idx=" << SPIdx <<endl; /////////////////////////////////////
				//Info <<"		Check: Ori=" << procSublayerProcOriCorr[thisProcNo][neighSublayerIdx] <<" ; Idx=" << procSublayerIdxCorr[thisProcNo][neighSublayerIdx] <<endl; /////////////////////////////////////
				procSublayerLinkCorr[thisProcNo][neighSublayerIdx] = infLinkWeight[i0][kS];
			}
		}
	}
}
Pstream::gatherList(procSublayerLinkCorr);
Pstream::scatterList(procSublayerLinkCorr);

//	// Extract values from neighbour processors
for(int pP=0; pP<procPatchIDs.size(); ++pP)
{
	int refp=procPatchRefNo[thisProcNo][pP];
	//Info << "Corr from "<< refp <<" ; Corr list size :" << procSublayerLinkCorr[refp].size() <<endl; /////////////////////////////////////
	if (procSublayerLinkCorr[refp].size()>0)
	{
		for (int nCorr=0; nCorr<procSublayerLinkCorr[refp].size(); ++nCorr)
		{
			//Info <<"	n corr: " << nCorr  << endl; /////////////////////////////////////
			//Info <<"	n corr " << nCorr << " to processor "<< procSublayerProcOriCorr[refp][nCorr] << endl; /////////////////////////////////////
			if (procSublayerProcOriCorr[refp][nCorr] == thisProcNo)
			{
				int idxCorr = procSublayerIdxCorr[refp][nCorr];
				checkInfLink[sublayerIDs[idxCorr]] += procSublayerLinkCorr[refp][nCorr];
			}
		}
	}
}

//Patch correction
List< List<scalar> > procPatchSaaConcInfArea(totProcNo);
List<scalar> procPatchSaaConcInfAreaP(nProcElements);
procPatchSaaConcInfArea[thisProcNo] = procPatchSaaConcInfAreaP;
Pstream::gatherList(procPatchSaaConcInfArea);
Pstream::scatterList(procPatchSaaConcInfArea);
