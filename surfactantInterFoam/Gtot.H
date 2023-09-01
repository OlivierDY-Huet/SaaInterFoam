{

//Info << "Gtot " << endl;


			
//Gradient in X, Y and Z directions
//		Internal values
scalar lim = 0.5;
forAll(alpha1().mesh().magSf(),facei)
{
	if (alpha1().mesh().isInternalFace(facei))
	{
		label celliO = alpha1().mesh().faceOwner()[facei];
		label celliN = alpha1().mesh().faceNeighbour()[facei];
		vector O2N(alpha1().mesh().C()[celliN]-alpha1().mesh().C()[celliO]);
		scalar distfi(mag(O2N));
		vector nfi(O2N/distfi);
			
		scalar Gfi = (psi[celliN]-psi[celliO])/distfi;
		
		if(nfi[0]>lim)
		{
			aa[celliO] = Gfi*nfi[0];
			bb[celliN] = Gfi*nfi[0];
		}
		else if(nfi[0]<-lim)
		{
			aa[celliN] = Gfi*nfi[0];
			bb[celliO] = Gfi*nfi[0];
		}
		
		if(nfi[1]>lim)
		{
			cc[celliO] = Gfi*nfi[1];
			dd[celliN] = Gfi*nfi[1];
		}
		else if(nfi[1]<-lim)
		{
			cc[celliN] = Gfi*nfi[1];
			dd[celliO] = Gfi*nfi[1];
		}
		
		if(nfi[2]>lim)
		{
			ee[celliO] = Gfi*nfi[2];
			ff[celliN] = Gfi*nfi[2];
		}
		else if(nfi[2]<-lim)
		{
			ee[celliN] = Gfi*nfi[2];
			ff[celliO] = Gfi*nfi[2];
		}
	}
	
}

aa.correctBoundaryConditions();
bb.correctBoundaryConditions();
cc.correctBoundaryConditions();
dd.correctBoundaryConditions();
ee.correctBoundaryConditions();
ff.correctBoundaryConditions();

//		Processor patch values
const fvPatchList& patches = mesh.boundary();
forAll(patches, patchi)
{
	const fvPatch& curPatch = patches[patchi];
	if (isType<processorFvPatch>(curPatch))
	{
		vectorField nf = curPatch.nf(); //
		
		const fvPatchScalarField& cellWidthp = cellWidth.boundaryField()[patchi];
		scalarField cellWidthpi = cellWidthp.patchInternalField();
		scalarField distfi(cellWidthpi/2.0+cellWidthp/2.0);
		
		const fvPatchScalarField& psip = psi.boundaryField()[patchi];
		scalarField psipi = psip.patchInternalField();
		
		scalarField Gf = (psip-psipi)/distfi;
		
		forAll(curPatch,facei)
		{	
			int celliO = curPatch.faceCells()[facei]; //Face patch is always owned by the cell next to the patch
			vector nfi = nf[facei];
			scalar Gfi = Gf[facei];
			
			if(nfi[0]>lim)
			{
				aa[celliO] = Gfi*nfi[0];
			}
			else if(nfi[0]<-lim)
			{
				bb[celliO] = Gfi*nfi[0];
			}
			
			if(nfi[1]>lim)
			{
				cc[celliO] = Gfi*nfi[1];
			}
			else if(nfi[1]<-lim)
			{
				dd[celliO] = Gfi*nfi[1];
			}
			
			if(nfi[2]>lim)
			{
				ee[celliO] = Gfi*nfi[2];
			}
			else if(nfi[2]<-lim)
			{
				ff[celliO] = Gfi*nfi[2];
			}
			
		
		}
	}
}

aa.correctBoundaryConditions();
bb.correctBoundaryConditions();
cc.correctBoundaryConditions();
dd.correctBoundaryConditions();
ee.correctBoundaryConditions();
ff.correctBoundaryConditions();



// Oposite signal

aaMod = aa;
bbMod = bb;
ccMod = cc;
ddMod = dd;
eeMod = ee;
ffMod = ff;

forAll(psi,celli)
{
	if(aa[celli]*bb[celli]<0.0)
	{
		scalar newG = aa[celli]+bb[celli];
		aaMod[celli] = newG;
		bbMod[celli] = newG;
	}
	if(cc[celli]*dd[celli]<0.0)
	{
		scalar newG = cc[celli]+dd[celli];
		ccMod[celli] = newG;
		ddMod[celli] = newG;
	}
	if(ee[celli]*ff[celli]<0.0)
	{
		scalar newG = ee[celli]+ff[celli];
		eeMod[celli] = newG;
		ffMod[celli] = newG;
	}
}



const volScalarField::Boundary& abf = alpha1.boundaryField();
forAll(patches,patchi)
{
	if (isA<alphaContactAngleTwoPhaseFvPatchScalarField>(abf[patchi]))
	{
		const fvPatchScalarField& aap = aa.boundaryField()[patchi];
		const fvPatchScalarField& bbp = bb.boundaryField()[patchi];
		const fvPatchScalarField& ccp = cc.boundaryField()[patchi];
		const fvPatchScalarField& ddp = dd.boundaryField()[patchi];
		const fvPatchScalarField& eep = ee.boundaryField()[patchi];
		const fvPatchScalarField& ffp = ff.boundaryField()[patchi];
		
		const fvPatch& curPatch = patches[patchi];
		vectorField nf = curPatch.nf();
		
		forAll(curPatch,facei)
		{
			int faceCelli = curPatch.faceCells()[facei];
			vector nfi = nf[facei];
			
			if(nfi[0]>lim)
			{ 	
				aaMod[faceCelli]=bb[faceCelli];
			}
			if(nfi[0]<-lim)
			{
				bbMod[faceCelli]=aa[faceCelli];
			}
		
			if(nfi[1]>lim)
			{
				ccMod[faceCelli]=dd[faceCelli];
			}
			if(nfi[1]<-lim)
			{
				ddMod[faceCelli]=cc[faceCelli];
			}
			
			if(nfi[2]>lim)
			{
				eeMod[faceCelli]=ff[faceCelli];
			}
			if(nfi[2]<-lim)
			{
				ffMod[faceCelli]=ee[faceCelli];
			}
		}
	}
}

aaMod.correctBoundaryConditions();
bbMod.correctBoundaryConditions();
ccMod.correctBoundaryConditions();
ddMod.correctBoundaryConditions();
eeMod.correctBoundaryConditions();
ffMod.correctBoundaryConditions();



//Gtot
/*
Gin = Foam::sqrt(
		max(pow(max(bbMod,double(0)),2),pow(min(aaMod,double(0)),2))
	  + max(pow(max(ddMod,double(0)),2),pow(min(ccMod,double(0)),2))
	  + max(pow(max(eeMod,double(0)),2),pow(min(ffMod,double(0)),2))
		);
Gout = Foam::sqrt(
		max(pow(min(bbMod,double(0)),2),pow(max(aaMod,double(0)),2))
	  + max(pow(min(ddMod,double(0)),2),pow(max(ccMod,double(0)),2))
	  + max(pow(min(eeMod,double(0)),2),pow(max(ffMod,double(0)),2))
		);
Gin.correctBoundaryConditions();
Gout.correctBoundaryConditions();


Gtot = pos0(alpha1-0.5) * Gin + pos(0.5-alpha1) * Gout;
Gtot.correctBoundaryConditions();
*/


forAll(alpha1,celli)
{
	//in
	if(pow(max(bbMod[celli],double(0)),2)>pow(min(aaMod[celli],double(0)),2))
	{
		Ginv[celli][0] = max(bbMod[celli],double(0));
	}
	else
	{
		Ginv[celli][0] = min(aaMod[celli],double(0));
	}
	if(pow(max(ddMod[celli],double(0)),2)>pow(min(ccMod[celli],double(0)),2))
	{
		Ginv[celli][1] = max(ddMod[celli],double(0));
	}
	else
	{
		Ginv[celli][1] = min(ccMod[celli],double(0));
	}
	if(pow(max(ffMod[celli],double(0)),2)>pow(min(eeMod[celli],double(0)),2))
	{
		Ginv[celli][2] = max(ffMod[celli],double(0));
	}
	else
	{
		Ginv[celli][2] = min(eeMod[celli],double(0));
	}
	
	//out
	if(pow(min(bbMod[celli],double(0)),2)>pow(max(aaMod[celli],double(0)),2))
	{
		Goutv[celli][0] = min(bbMod[celli],double(0));
	}
	else
	{
		Goutv[celli][0] = max(aaMod[celli],double(0));
	}
	if(pow(min(ddMod[celli],double(0)),2)>pow(max(ccMod[celli],double(0)),2))
	{
		Goutv[celli][1] = min(ddMod[celli],double(0));
	}
	else
	{
		Goutv[celli][1] = max(ccMod[celli],double(0));
	}
	if(pow(min(ffMod[celli],double(0)),2)>pow(max(eeMod[celli],double(0)),2))
	{
		Goutv[celli][2] = min(ffMod[celli],double(0));
	}
	else
	{
		Goutv[celli][2] = max(eeMod[celli],double(0));
	}
}
Ginv.correctBoundaryConditions();
Goutv.correctBoundaryConditions();

Gtotv = pos0(alpha1-0.5) * Ginv + pos(0.5-alpha1) * Goutv;
Gtotv.correctBoundaryConditions();

}
