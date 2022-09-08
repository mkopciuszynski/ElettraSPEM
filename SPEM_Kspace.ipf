#pragma rtGlobals=3     // Use modern global access method.
// This file contains procedures used to transform the data to K-space


// This function performs the transformation to K space of a single 2D map
// It is not accurate, the k values are accurate only for scans done around Gamma points (P~90)
function Transform2Dscan2K(srcName)
	String srcName

	wave data2D = $srcName
	/////////////////////////
	Variable i,j,k
	Variable numOfE,numOfNu   // number of Energy/detector angles points
	/////////////////////////
	// get the goniometer data
	SetDataFolder ksDataFolder
	Wave/T txtAttr=$ksTextAttributes
	Variable PP,TT
	PP=str2num(txtAttr[%$srcName][%$ksaP])
	TT=str2num(txtAttr[%$srcName][%$ksaT])
	/////////////////////////////////////////////
	Variable startE, deltaE, deltaNu, startNu
	deltaE = DimDelta($srcName,0)
	startE = DimOffset($srcName,0)
	numOfE = DimSize($srcName,0)
	numOfNu= DimSize($srcName,1)
	deltaNu= DimDelta($srcName,1)
	startNu= DimOffset($srcName,1)
	/////////////////////////////////////////////
	nvar gEf27 = root:SPMData:gEf27
	nvar gEf74 = root:SPMData:gEf74
	Variable Ef = gEf27
	if (startE>30)
		Ef = gEf74
	endif

	/////////////////////////////////////////////
	String dataKname = srcName +"_K"
	make /O/N=(DimSize(data2D,0),DimSize(data2D,1)) $dataKName
	wave dataK = $dataKName
	/////////////////////////////////////////////
	variable Kx1=Ximp(startE,PP,TT,startNu,0,0)
	variable Ky1=Yimp(startE,PP,TT,startNu,0,0)
	variable Kx2=Ximp(startE,PP,TT,startNu+numOfNu*deltaNu,0,0)
	variable Ky2=Yimp(startE,PP,TT,startNu+numOfNu*deltaNu,0,0)
	/////////////////////////////////////////////
	variable kxstart,kystart,kxend,kyend
	kxstart = min(kx1,kx2)
	kxend = max(kx1,kx2)

	kyend = max(ky1,ky2)
	kystart = min(ky1,ky2)
	/////////////////////////////////////////////
	variable Ek
	Variable kx,ky
	kx = (kxstart+kxend)*0.5
	if (kxend-kxstart>0.1)
		print("Transformation to K is not acurate!!!")
	endif
	SetScale/I y, Kystart, Kyend, dataK
	SetScale/P x, startE,deltaE, dataK
	/////////////////////////////////////////////
	for (k=0; k<numOfE;k+=1)
		Ek = startE+deltaE*k
		for(j=0;j<DimSize(dataK,1);j+=1)
			ky = Kystart+DimDelta(dataK,1)*j

			dataK[k][j]=Interp2D(data2D,Ek,nufromKXKYv2(kx,ky,TT,Ek))
		endfor
	endfor
	SetScale/P x, startE-Ef,deltaE, dataK

	//SetScale/P y, Dimoffset(sliceEc_K,1),Dimdelta(sliceEc_K,1), "Kx", dataK
	//SetScale/P y, Dimoffset(sliceEc_K,1),Dimdelta(sliceEc_K,1), "Ky", dataK
	//SetScale/P x, startE-Ef, deltaE, "Ek - Ef", dataK
	plotEk(dataKname)
end function



// TODO Recognise if the merging is done in reverese order
// This function merges two maps (the bottom could be empty)
function Merge(bottomDataName,topDataName)
	String bottomDataName,topDataName

	// TODO
	//	if (WaveExists($topDataName) && WaveExists($bottomDataName))

	Wave topData = $topDataName
	Wave bottomData = $bottomDataName

	Variable i,j,k,l

	Variable kxMin,kxMax,kyMin,kyMax,kxStep,kyStep
	Variable tkx,tky,bxi,byi
	kxMin = DimOffset(bottomData,0)
	kxStep = DimDelta(bottomData,0)
	kxMax = DimOffset(bottomData,0)+DimSize(bottomData,0)*kxStep
	kyMin = DimOffset(bottomData,1)
	kyStep = DimDelta(bottomData,1)
	kyMax = DimOffset(bottomData,1)+DimSize(bottomData,1)*kyStep

	// Project the top layer onto the same kx ky coordinates as the bottom layer
	Duplicate/O bottomData topProData
	topProData[][][] = NaN

	for(k=0; k<DimSize(topProData,2);k+=1)          // energy
		for(i=0;i<DimSize(topProData,0); i+=1)			// kx
			for(j=0;j<DimSize(topProData,1);j+=1)		// ky
				tkx = Kxmin+i*Kxstep
				tky = Kymin+j*Kystep
				if     ((dimoffset(topdata,0)<tkx) && (dimoffset(topdata,0)+dimsize(topdata,0)*dimdelta(topdata,0)-0.01)>tkx)
					if ((dimoffset(topdata,1)<tky) && (dimoffset(topdata,1)+dimsize(topdata,1)*dimdelta(topdata,1)-0.01)>tky)
						topProData[i][j][k] = topData(tkx)(tky)[k]
					endif
				endif
			endfor
		endfor
	endfor

	Make/O/N=(DimSize(bottomData,1)) valBottom
	Make/O/N=(DimSize(bottomData,1)) valTop
	Make/B/O/N=(DimSize(bottomData,1)) maskTop
	Make/B/O/N=(DimSize(bottomData,1)) maskBottom
	Make/B/O/N=(DimSize(bottomData,1)) maskOverlap
	Make/B/O/N=(DimSize(bottomData,1)) maskSum

	variable numOverlap,delta,indStart

	for(k=0; k<DimSize(bottomData,2);k+=1)			// energy
		for(i=0;i<DimSize(bottomData,0); i+=1)		// kx
			// get one collumn of the bottom and top data
			valTop[]    = topProData[i][p][k]
			valBottom[] = bottomData[i][p][k]
			// data present means 1
			MatrixOp/O maskTop	 = -0.5*(numType(valTop)-2)
			MatrixOp/O maskBottom = -0.5*(numType(valBottom)-2)
			MatrixOp/O maskSum = (maskTop || maskBottom)
			MatrixOp/O maskOverlap = (maskTop)*(MaskBottom)
			numOverlap = sum(maskOverlap) // How many points are overlaping
			delta = 1/numOverlap
			// Find the starting point of the overlaping
			FindLevel/Q maskOverlap 1
			indStart = V_LevelX
			// Create linear decay on both datasets in the overlaping region
			for(j=0;j<DimSize(bottomData,1);j+=1)   // ky
				if (maskOverlap[j]>0)
					valTop[j]    *= delta*(j-indStart)
					valBottom[j] *= 1-delta*(j-indStart)
				endif
			endfor

			MatriXOp/O valBottom = replaceNaNs(valBottom,0)
			MatriXOp/O valTop = replaceNaNs(valTop,0)
			MatrixOp/O valBottom = valBottom + valTop
			MatriXOp/O valBottom = setNaNs(valBottom,maskSum-1)

			bottomData[i][][k] = valBottom[q]
		endfor
	endfor

	//KillWaves valTop,valBottom,maskTop,maskBottom,maskOverlap,maskSum,topProData
end function


/////////////////////////////////////////////
// P - goniometer angle perpendicular to the slit -- for P90 the 90 deg means normal emission
// T - goniometer angle parallel to the slit
// Nu - angles on the detector about +/- 8 deg

// In the area measurements T is constant, P is changed

// T0 - value for normal emission
// P90 - value for normal emission
/////////////////////////////////////////////

//////////////////////////////////////////////////
// Transform and merge the whole region
// This function recognizes the type of measurements: region or line
// The user can chose if the data should be merged
function TransformMerge2K(srcName)
	String srcName
	String tName,newName,stemp
	/////////////////////////////////////////////
	// Get the names of all waves from this 3D measurement
	String st,list              // temp strig / list of data
	tName = "SMPM"+srcName[3,strSearch(srcName,"T",0)-1]
	list = WaveList(tName+"*",";","")
	Wave/T listWv = ListToTextWave(list,";")    // Convert to text wave
	Variable listLen = DimSize(listWv,0)
	/////////////////////////
	// get the goniometer data
	Make/O/N=(2,listLen) goniometerPT   // prepare 2 rows wave for phi and theta values
	String fldrSav=GetDataFolder(1)
	SetDataFolder ksDataFolder
	Wave/T txtAttr=$ksTextAttributes
	/////////////////////////
	Variable i,j,k,l
	Variable numOfAreas=1                   // Number of measured areas
	Variable numOfP=1                       // how many P points in one area
	Variable numOfE,numOfNu                 // number of Energy/detector angles points
	/////////////////////////
	for(i=0;i<listLen;i+=1)
		st=listWv[i]
		goniometerPT[0][i]=str2num(txtAttr[%$st][%$ksaP])
		goniometerPT[1][i]=str2num(txtAttr[%$st][%$ksaT])
		if (i>=1) // check how many areas were measured
			if (abs(goniometerPT[1][i]-goniometerPT[1][i-1])>1)     // detect if theta steps are bigger that 1 degree
				numOfAreas+=1
			endif
			if(numOfAreas==1)                                               // count the number of points in a single area
				numOfP+=1
			endif
		endif
	endfor
	/////////////////////////////////////////////
	nvar gEf27 = root:SPMData:gEf27
	nvar gEf74 = root:SPMData:gEf74
	nvar gT0 = root:SPMData:gT0
	nvar gP90 = root:SPMData:gP90
	Variable Ef = gEf27
	if ((DimOffset($srcName,2)+Dimsize($srcName,2)*Dimdelta($srcName,2))>60)
		Ef=gEf74
	endif

	Variable P90=gP90,T0=gT0, dgonP=0
	String transType,mergeAreas
	Prompt transType ,"Transformation type",popup,"simple;complete"
	Prompt mergeAreas,"Merge areas?",popup,"Yes;NO"
	Prompt Ef       ,"Ef"
	Prompt P90      ,"P for normal emission"
	Prompt T0       ,"T for normal emission"
	Prompt dgonP    ,"Goniometer misalignment in P to add"
	if (numOfAreas>1)
		Doprompt    "Measurements parameters",transType,mergeAreas,Ef,P90,T0,dgonP
	else
		Doprompt    "Measurements parameters",transType,Ef,P90,T0,dgonP
		mergeAreas = "NO"
	endif
	gT0 = T0
	gP90 = P90

	// User not canceled
	if (V_flag==0)

		/////////////////////////////////////////////
		// Progress panel
		/////////////////////////////////////////////
		NewPanel /N=ProgressPanel /W=(300,200,800,260)
		TitleBox   winmsg,  pos={20,10},   size={170,15},title="", frame=0
		ValDisplay valdisp0,pos={20,30},size={400,20},limits={0,(numOfP),0},barmisc={0,0}
		ValDisplay valdisp0,value= _NUM:0,mode=3,highColor=(0,65535,0)
		Button bStop,pos={440,30},size={50,20},title="Stop"
		DoUpdate /W=ProgressPanel /E=1
		/////////////////////////////////////////////
		/////////////////////////////////////////////
		// in Merge == NO remove some data from the 3D scan
		// find the start and end of selected data
		if (numOfAreas>1)
			if (stringmatch(mergeAreas,"NO"))
				variable indNext,indStart,indNum
				// T for measured part of the region
				variable srcT = str2num(srcName[strSearch(srcName,"T",0)+1,20])
				for (i=0;i<DimSize(goniometerPT,1);i+=1)
					// find area in region according to the T value
					if (abs(goniometerPT[1][i]-srcT)<0.1)
						if (indNext==0)
							indStart = i; indNext = 1
						else // count the points
							indNum+=1
						endif
					endif
				endfor
				DeletePoints /M=1 0,indStart,goniometerPT
				DeletePoints /M=1 indNum+1,dimSize(goniometerPT,1),goniometerPT
				DeletePoints /M=0 0,indStart,listWV
				DeletePoints /M=0 indNum+1,dimSize(listWv,0),listWv
				numOfAreas = 1
			endif

		endif
		/////////////////////////////////////////////
		if ((numOfAreas*numOfP)!=dimSize(goniometerPT,1))
			print "Missing data files"
			return 0
		endif
		listLen = DimSize(listWv,0)
		/////////////////////////////////////////////
		Make/O/N=(DimSize(goniometerPT,1)) goniometerT
		goniometerT[]=goniometerPT[1][p]
		Variable deltaP, startP, startE, deltaE, deltaNu, startNu
		tName = listWv[0]
		startP = goniometerPT[0][0]
		deltaP = goniometerPT[0][1]-startP
		deltaE = DimDelta($tName,0)
		startE = DimOffset($tName,0)
		numOfE = DimSize($tName,0)
		numOfNu=DimSize($tName,1)
		deltaNu=DimDelta($tName,1)
		startNu=DimOffset($tName,1)
		// Prepare 4D matrix: P/Nu/Ek/regions in Areas
		Make/O/N=(numOfP,numOfNu,numOfE,numOfAreas) Ang4D
		SetScale/P x, startP, deltaP, "P (deg)" , Ang4D
		SetScale/P y, startNu,deltaNu,"nu (deg)", Ang4D
		SetScale/P z, startE, deltaE,"Ek (eV)",   Ang4D
		Ang4D[][][][]=NaN
		// load the data directly from 2D waves NOT FROM 3D angle maps!
		print "Data loaded directly form 2D waves"
		for (l=0;l<numOfAreas;l+=1)
			for (i=0;i<numOfP;i+=1)
				if (strlen(listWv[i])>0)
					tName = listWv[i+l*numOfP]
					Wave data2D = $tName
					Ang4D[i][][][l] = data2D[r][q]
				endif
			endfor
		endfor
		/////////////////////////////////////////////
		make/O/N=(numOfAreas) Kxmi
		make/O/N=(numOfAreas) Kxma
		make/O/N=(numOfAreas) Kymi
		make/O/N=(numOfAreas) Kyma
		Variable Kxmin, Kxmax, Kymin, Kymax, Kxstep, Kystep
		variable Ek

		TitleBox   winmsg,title="Transform to K space",win=ProgressPanel
		ValDisplay valdisp0,limits={0,(numOfAreas*numOfE),0},win=ProgressPanel

		for (l=0; l<numOfAreas; l+=1)
			Make/O/N=(numOfP) goniometerT
			goniometerT[] = goniometerPT[1][p+l*numOfP]
			srcT = mean(goniometerT)
			// Make a 2D wave for Econst slice
			Make/O/N=(Dimsize(Ang4D,0),Dimsize(Ang4D,1)) sliceEc_A
			SetScale/P x, Dimoffset(Ang4D,0)+dgonP, Dimdelta(Ang4D,0), sliceEc_A
			SetScale/P y, Dimoffset(Ang4D,1), Dimdelta(Ang4D,1), sliceEc_A

			Variable kx1,kx2,kx3,kx4,ky1,ky2,ky3,ky4,kxstart,kxend,kystart,kyend
			Variable T=goniometerT[0]
			Kx1=Ximp(Ef,startP,T,startNu,T0,P90-90)
			Ky1=Yimp(Ef,startP,T,startNu,T0,P90-90)
			Kx2=Ximp(Ef,startP,T,startNu+numOfNu*deltaNu,T0,P90-90)
			Ky2=Yimp(Ef,startP,T,startNu+numOfNu*deltaNu,T0,P90-90)
			T = goniometerT[numOfP-1]
			Kx3=Ximp(Ef,startP+numOfP*deltaP,T,startNu,T0,P90-90)
			Ky3=Yimp(Ef,startP+numOfP*deltaP,T,startNu,T0,P90-90)
			Kx4=Ximp(Ef,startP+numOfP*deltaP,T,startNu+numOfNu*deltaNu,T0,P90-90)
			Ky4=Yimp(Ef,startP+numOfP*deltaP,T,startNu+numOfNu*deltaNu,T0,P90-90)

			Kxstart = min(Kx1,Kx2,Kx3,Kx4)
			Kxend   = max(Kx1,Kx2,Kx3,Kx4)
			Kystart = min(Ky1,Ky2,Ky3,Ky4)
			Kyend   = max(Ky1,Ky2,Ky3,Ky4)
			variable KYsize=numOfNu*(1+((Kyend-Kystart)/abs(Ky2-Ky1)))

			Make/O/N=(numOfP,KYsize)  sliceEc_K
			Make/O/N=(numOfP,KYsize,numOfE) Ang3D_K

			for (i=0; i<numOfE;i+=1)
				sliceEc_A[][] = Ang4D[p][q][i][l]
				Ek = startE+deltaE*i
				TransformEcAng2K("sliceEc_A","sliceEc_K",Ek,"goniometerT",T0,P90, Ef,transType)
				//////////////////////////////////////////////
				Ang3D_K[][][i]=sliceEc_K[p][q]
				//////////////////////////////////////////////
				ValDisplay valdisp0,value= _NUM:(i+1+l*numOfE),win=ProgressPanel
				DoUpdate /W=ProgressPanel
				if (V_Flag == 2)
					break  // Break if stop button was pressed

			endif
		endfor

		SetScale/P x, Dimoffset(sliceEc_K,0),Dimdelta(sliceEc_K,0), "Kx", Ang3D_K
		SetScale/P y, Dimoffset(sliceEc_K,1),Dimdelta(sliceEc_K,1), "Ky", Ang3D_K
		SetScale/P z, startE-Ef, deltaE, "Ek", Ang3D_K

		Kxmi[L]=Dimoffset(sliceEc_K,0)
		Kymi[L]=Dimoffset(sliceEc_K,1)
		Kxma[L]=Kxmi[L]+Dimdelta(sliceEc_K,0)*Dimsize(sliceEc_K,0)
		Kyma[L]=Kymi[L]+Dimdelta(sliceEc_K,1)*Dimsize(sliceEc_K,1)

		if (Kymi[L]<=Kymin)
			Kxstep=Dimdelta(sliceEc_K,0)
			Kystep=Dimdelta(sliceEc_K,1)
			Kymin=Kymi[L]
		endif
		if (numOfAreas>1)
			NewName="D3_"+num2str(L)+"_K"
		else
			NewName=srcName+"_K"
		endif
		duplicate/O Ang3D_K $NewName
	endfor

	TitleBox   winmsg,title="Merge in progress",win=ProgressPanel
	ValDisplay valdisp0,value= _NUM:0,limits={0,(numOfAreas),0},win=ProgressPanel
	DoUpdate /W=ProgressPanel



	if (stringmatch(mergeAreas,"Yes"))
		print "Merge the whole region"

		newName=srcName[0,strSearch(srcName,"T",0)-1]+"merged_K"

		Kxmin=wavemin(Kxmi)
		Kxmax=wavemax(Kxma)
		Kymax=wavemax(Kyma)
		Kymin=wavemin(Kymi)

		Make/O/N=(((Kxmax-Kxmin)/Kxstep), ((Kymax-Kymin)/Kystep), Dimsize(Ang3D_K, 2)) $newName

		SetScale/P x, Kxmin, Kxstep, "Kx", $newName
		SetScale/P y, Kymin, Kystep, "Ky", $newName
		SetScale/P z, DimOffset(Ang3D_K, 2), DimDelta(Ang3D_K,2), "Ek", $newName

		wave D3_merged_K = $newName
		D3_merged_K = NaN

		for(l=0;l<numOfAreas;l+=1)
			tName="D3_"+num2str(L)+"_K"

			Merge(newName,tName)

			ValDisplay valdisp0,value= _NUM:(l+1),win=ProgressPanel
			DoUpdate /W=ProgressPanel
			if( V_Flag == 2)
				break           // Break if stop button was pressed
		endif

	endfor
else


endif

DisplayEcMap(newName)
KillWindow ProgressPanel

/// TODO remove all data used only in this function
//KillWaves Kxmi,Kxma,Kymi,Kyma,goniometerT,goniometerPT,Ang3D_K,Ang4D
	endif
end function




/////////////////////////////////////////////////////////
// Transform single slice of constant energy angle map to K space
// ImageAngle/K - name of the input and output waves
// E - kinetic energy
// gonio - name of the wave containing goniometer T values for each column
// T0 / P90 - the gamma point
// Ef - the Fermi level
// transType = simple/complete
function TransformEcAng2K(ImageAngle, ImageK,E,gonio,T0,P90,Ef,transType)
	string ImageAngle, ImageK, gonio,transType
	variable E, T0, P90, Ef
	variable T
	Wave dataK=$ImageK, dataAn=$ImageAngle, gonioT=$gonio

	variable sizeP  =DimSize(dataAn,0)
	variable sizeNu =DimSize(dataAn,1)
	variable startP =Dimoffset(dataAn,0)
	variable startNu=Dimoffset(dataAn,1)
	variable stepP  =Dimdelta(dataAn,0)
	variable stepNu =Dimdelta(dataAn,1)

	T=gonioT[0]
	variable Kx1=Ximp(Ef,startP,T,startNu,T0,P90-90)
	variable Ky1=Yimp(Ef,startP,T,startNu,T0,P90-90)
	variable Kx2=Ximp(Ef,startP,T,startNu+sizeNu*stepNu,T0,P90-90)
	variable Ky2=Yimp(Ef,startP,T,startNu+sizeNu*stepNu,T0,P90-90)
	T=gonioT[sizeP-1]
	variable Kx3=Ximp(Ef,startP+sizeP*stepP,T,startNu,T0,P90-90)
	variable Ky3=Yimp(Ef,startP+sizeP*stepP,T,startNu,T0,P90-90)
	variable Kx4=Ximp(Ef,startP+sizeP*stepP,T,startNu+sizeNu*stepNu,T0,P90-90)
	variable Ky4=Yimp(Ef,startP+sizeP*stepP,T,startNu+sizeNu*stepNu,T0,P90-90)
	variable kxstart,kystart,kxend,kyend
	kxstart = min(kx1,kx2,kx3,kx4)
	kystart = min(ky1,ky2,ky3,ky4)
	kxend = max(kx1,kx2,kx3,kx4)
	kyend = max(ky1,ky2,ky3,ky4)

	SetScale/I x, Kxstart, Kxend, dataK
	SetScale/I y, Kystart, Kyend, dataK
	Variable kx,ky
	variable i, j
	for(i=0;i<DimSize(dataK,0);i+=1)
		for(j=0;j<DimSize(dataK,1);j+=1)
			kx = Kxstart+DimDelta(dataK,0)*i
			ky = Kystart+DimDelta(dataK,1)*j
			T=gonioT[sizeP-1-i]


			if (stringmatch(transType,"complete"))
				// complete transformation
				Wave nuP=Pandnufromk(E,T,T0,P90-90,kx,ky)
				dataK[i][j]=Interp2D(dataAn,nuP[1],nuP[0])
			else
				// simple tranformation
				T=gonioT[sizeP-1-i]-T0
				dataK[i][j]=Interp2D(dataAn,(PfromKXKY(kx,ky,T, E)),(nufromKXKY(kx,ky,T, E)))
			endif
		endfor
	endfor
end function


function/Wave Pandnufromk(E,T,offT,offP,kx,ky)
	variable E,T,offT,offP,kx,ky
	variable ind=0.000001, indi
	variable i=0
	variable nu,P, kxi, kyi, nu0, P0
	Make/O nuP={0,0}

	nu = nufromKXKY(kx,ky,T-offT,E)
	P  = PfromKXKY(kx,ky,T-offT,E)+offP
	nuP[0]=nu
	nuP[1]=P

	kxi=Ximp(E,P,T,nu,offT,offP)
	kyi=Yimp(E,P,T,nu,offT,offP)

	if (offP==0)
		return nuP
	else
		do
			nu0=nu
			P0=P
			nu=nunext(E,P0,T,nu0,offT,offP,kx,ky)
			P=Pnext(E,P0,T,nu0,offT,offP,kx,ky)

			kxi=Ximp(E,P,T,nu,offT,offP)
			kyi=Yimp(E,P,T,nu,offT,offP)

			indi=(kx-kxi)*(kx-kxi)+(ky-kyi)*(ky-kyi)
			i+=1

		while ((indi>ind)&&(i<5))
		nuP[0]=nu
		nuP[1]=P
		return nuP
	endif

end



function TransformEcAng2K_new(ImageAngle, ImageK,E,gonio,T0,P90,Ef)
	string ImageAngle, ImageK, gonio
	variable E, T0, P90, Ef
	variable T
	Wave dataK=$ImageK, dataAn=$ImageAngle, gonioT=$gonio
	variable sizeP  =DimSize(dataAn,0)
	variable sizeNu =DimSize(dataAn,1)
	variable startP =Dimoffset(dataAn,0)
	variable startNu=Dimoffset(dataAn,1)
	variable stepP  =Dimdelta(dataAn,0)
	variable stepNu =Dimdelta(dataAn,1)

	Make/O/N=(sizeP*sizeNu,3) sampleTriplet

	Variable i,j,tP,tNu
	for(i=0;i<sizeP;i+=1)
		for(j=0;j<sizeNu;j+=1)
			T=gonioT[sizeP-1-i]
			tP = startP+stepP*i
			tNu = startNu + stepNu*j
			sampleTriplet[j+i*sizeNu][0] = xImp(E,tP,T,tNu,0,0)
			sampleTriplet[j+i*sizeNu][1] = yImp(E,tP,T,tNu,0,0)
			//sampleTriplet[j+i*sizeNu][0] = kxFromPTnu(E,tP,T,tNu,0,0)
			//sampleTriplet[j+i*sizeNu][1] = kyFromPTnu(E,tP,T,tNu,0,0)
			sampleTriplet[j+i*sizeNu][2] = dataAn[i][j]
		endfor
	endfor

	ImageInterpolate/RESL={300,300}/DEST=firstImage voronoi, sampleTriplet


	//ImageInterpolate/S={1.5,0.01,100,1.5,0.01,100} /K={1,0, 1, 1} /DEST=firstImage Kriging, sampleTriplet


	//Make /O /N=(100,100) dataMat=0
	//SetScale x,0,2,dataMat
	//SetScale y,0,2,dataMat
	//Duplicate /O dataMat,countMat
	//ImageFromXYZ /AS sampleTriplet, dataMat, countMat

end function


// old functions

function Ximp(E,P,T,nu,offT,offP)
	variable E, T, offT, nu, P, offP
	P*=Pi/180;T*=Pi/180;nu*=Pi/180;offT*=Pi/180;offP*=Pi/180
	return 0.512*sqrt(E)*(cos(nu)*(Cos(P)*Cos(offP)+Sin(P)*Sin(offP)*Cos(T-offT))-Sin(nu)*Sin(offP)*sin(T-offT))
end

function Yimp(E,P,T,nu,offT,offP)
	variable E, T, offT, nu, P, offP
	P*=Pi/180;T*=Pi/180;nu*=Pi/180;offT*=Pi/180;offP*=Pi/180
	return 0.512*sqrt(E)*(Cos(nu)*Sin(P)*Sin(T-offT)+Sin(nu)*Cos(T-offT))
end
function nufromKXKYv2(Ximp,Yimp,T,E)
	variable Ximp, Yimp, T, E
	Variable nu1, nu2
	Ximp /=0.512*sqrt(E)
	Yimp /=0.512*sqrt(E)
	if ((Ximp*Ximp+Yimp*Yimp)>1)
		return NaN
	endif
	if (T==0)
		return 180/pi*asin(Yimp)
	else
		nu1=Yimp*Cos(pi/180*T)+Sin(pi/180*T)*sqrt(1-Ximp*Ximp-Yimp*Yimp)
		nu2=Yimp*Cos(pi/180*T)-Sin(pi/180*T)*sqrt(1-Ximp*Ximp-Yimp*Yimp)
		nu1=180/pi*asin(nu1)
		nu2=180/pi*asin(nu2)
		Variable Ky1, Ky2
		Variable P
		P=180/pi*acos(Ximp/Cos(pi/180*nu1))
		Ky1=Yimpulse(E,P,T,nu1)
		P=180/pi*acos(Ximp/Cos(pi/180*nu2))
		Ky2=Yimpulse(E,P,T,nu2)
		Ky1/=0.512*sqrt(E)
		Ky2/=0.512*sqrt(E)
		return nu2
	endif
end

function nufromKXKY(Ximp,Yimp,T,E)
	variable Ximp, Yimp, T, E
	Variable nu1, nu2
	Ximp /=0.512*sqrt(E)
	Yimp /=0.512*sqrt(E)
	if ((Ximp*Ximp+Yimp*Yimp)>1)
		return NaN
	endif
	if (T==0)
		return 180/pi*asin(Yimp)
	else
		nu1=Yimp*Cos(pi/180*T)+Sin(pi/180*T)*sqrt(1-Ximp*Ximp-Yimp*Yimp)
		nu2=Yimp*Cos(pi/180*T)-Sin(pi/180*T)*sqrt(1-Ximp*Ximp-Yimp*Yimp)
		nu1=180/pi*asin(nu1)
		nu2=180/pi*asin(nu2)
		Variable Ky1, Ky2
		Variable P
		P=180/pi*acos(Ximp/Cos(pi/180*nu1))
		Ky1=Yimpulse(E,P,T,nu1)
		P=180/pi*acos(Ximp/Cos(pi/180*nu2))
		Ky2=Yimpulse(E,P,T,nu2)
		Ky1/=0.512*sqrt(E)
		Ky2/=0.512*sqrt(E)
		if (abs(Ky1-Yimp)<0.001)
			return nu1
		elseif (abs(Ky2-Yimp)<0.001)
			return nu2
		else
			return NaN
		endif
	endif
end

function PfromKXKY(Ximpulse,Yimpulse,T, E)
	variable Ximpulse, Yimpulse, T, E
	variable nu
	nu=nufromKXKY(Ximpulse,Yimpulse,T, E)
	Ximpulse/=0.512*sqrt(E)
	if (numtype(nu)==0)
		return 180/pi*acos(Ximpulse/Cos(pi/180*nu))
	else
		return NaN
	endif
end


function nunext(E,P,T,nu,offT,offP,kx,ky)
	variable E, T, offT, nu, P, offP, kx,ky
	variable nunxt, kxi, kyi, determinant

	kxi=Ximp(E,P,T,nu,offT,offP)
	kyi=Yimp(E,P,T,nu,offT,offP)
	determinant=(dfxdnu(E,P,T,nu,offT,offP)*dfydP(E,P,T,nu,offT,offP))-(dfxdP(E,P,T,nu,offT,offP)*dfydnu(E,P,T,nu,offT,offP))
	nunxt=(Pi/180)*nu+(((kx-kxi)*dfydP(E,P,T,nu,offT,offP)-(ky-kyi)*dfxdP(E,P,T,nu,offT,offP))/determinant)

	return 180*nunxt/Pi
end


function Pnext(E,P,T,nu,offT,offP,kx,ky)
	variable E, T, offT, nu, P, offP, kx,ky
	variable Pnxt, kxi, kyi, determinant

	kxi=Ximp(E,P,T,nu,offT,offP)
	kyi=Yimp(E,P,T,nu,offT,offP)

	determinant=(dfxdnu(E,P,T,nu,offT,offP)*dfydP(E,P,T,nu,offT,offP))-(dfxdP(E,P,T,nu,offT,offP)*dfydnu(E,P,T,nu,offT,offP))

	Pnxt=(Pi/180)*P+(((ky-kyi)*dfxdnu(E,P,T,nu,offT,offP)-(kx-kxi)*dfydnu(E,P,T,nu,offT,offP))/determinant)
	return 180*Pnxt/Pi
end


Function Yimpulse(E, P, T,nu) //Ky from angles
	Variable E,P,T,nu
	Variable Ky
	Ky=0.512*Sqrt(E)*(Sin(pi/180*nu)*Cos(pi/180*T)+Sin(pi/180*P)*Sin(pi/180*T)*Cos(pi/180*nu))
	return Ky
End

function dfxdnu(E,P,T,nu,offT,offP)
	variable E, T, offT, nu, P, offP
	P*=Pi/180; T*=Pi/180; nu*=Pi/180; offT*=Pi/180; offP*=Pi/180;
	return -0.512*sqrt(E)*(Sin(nu)*(Cos(P)*Cos(offP)+Sin(P)*Sin(offP)*Cos(T-offT))+Cos(nu)*Sin(offP)*sin(T-offT))
end

function dfydnu(E,P,T,nu,offT,offP)
	variable E, T, offT, nu, P, offP
	P*=Pi/180;  T*=Pi/180;  nu*=Pi/180; offT*=Pi/180; offP*=Pi/180
	return 0.512*sqrt(E)*(Cos(nu)*Cos(T-offT)-Sin(nu)*Sin(P)*Sin(T-offT))
end

function dfydP(E,P,T,nu,offT,offP)
	variable E, T, offT, nu, P, offP
	P*=Pi/180;  T*=Pi/180;  nu*=Pi/180; offT*=Pi/180; offP*=Pi/180
	return 0.512*sqrt(E)*Cos(nu)*Cos(P)*Sin(T-offT)
end

function dfxdP(E,P,T,nu,offT,offP)
	variable E, T, offT, nu, P, offP
	P*=Pi/180;  T*=Pi/180;  nu*=Pi/180; offT*=Pi/180; offP*=Pi/180
	return 0.512*sqrt(E)*(cos(nu)*(Cos(P)*Sin(offP)*Cos(T-offT)-Sin(P)*Cos(offP)))
end


