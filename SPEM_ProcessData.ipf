#pragma rtGlobals=3   // Use modern global access method.

// Calculate the 2nd derivative of the E(k) map
Function Calc2ndDiff(dataName)
	String dataName

	wave srcData = $dataName
	//String dstDataName = dataName + "_DIF"
	Smooth 30, srcData
	srcData[][] = -srcData[p][q]
	Differentiate/DIM=0 /EP=1 srcData /D=srcData
	Differentiate/Dim=0 /EP=1 srcData /D=srcData
	DeletePoints 0,5, srcData
	DeletePoints DimSize(srcData,0)-5,5, srcData
	//DeletePoints 0,4, srcData

End


// Create the new dataset for spectro image integrated over the energy/energy+angle range
Function/T ImageSum(dataName)
	String dataName
	String dataSumName=dataName+"_sum"

	if(WaveDims($dataName)==4) // 4D image sum over the first and the second dimension
		String dataSumE = dataName+"_sumE"
		ReduceDim(dataName,0,dataSumE)
		dataName = dataSumE
		ReduceDim(dataName,0,dataSumName)
		killwaves $dataSumE
	endif
	if(WaveDims($dataName)==3)	// sum data over the first dimension
		ReduceDim(dataName,0,dataSumName)
	endif

	return dataSumName
end function


// Reduce the data dimension by summation over the selected dimension (Dim)
// Start, Finish define the range of summation
Function ReduceDim(dataName,Dim,dataSumName,[Start,Finish])
	String dataName,dataSumName
	Variable Dim,Start,Finish
	Variable i,j,k,m
	Variable wholeRange
	// Check if optional parameters were provided, otherwise the whole range will be used
	if( ParamIsDefault(Start) && ParamIsDefault(Finish) )
		wholeRange = 1
	endif

	Variable dataDim=WaveDims($dataName)
	Make/O dimList = {0,1,2,3}

	if(WaveExists($dataName) && strlen(dataSumName) && (dataDim>1) && (Dim<dataDim))
		WAVE srcData=$dataName
		// Make a copy of src Data for dimension trimming
		Duplicate /O srcData, srcCopy
		switch (dataDim)// data dimension
			case 2:
				DeletePoints dim,1, dimList
				if (!wholeRange)
					if (dim==0)
						duplicate /O /RMD=[start,finish][][] srcData,srcCopy
					elseif (dim==1)
						duplicate /O /RMD=[][start,finish][] srcData,srcCopy
					endif
				endif
				SumDimension /D=(dim)/DEST=$dataSumName srcCopy
				CopyWaveAttributes(dataName,dimList[0],dataSumName,0)
				break
			case 3:
				DeletePoints dim,1, dimList
				if (!wholeRange)
					if (dim==0)
						duplicate /O /RMD=[start,finish][][] srcData,srcCopy
					elseif (dim==1)
						duplicate /O /RMD=[][start,finish][] srcData,srcCopy
					elseif (dim==2)
						duplicate /O /RMD=[][][start,finish] srcData,srcCopy
					endif
				endif
				SumDimension /D=(dim)/DEST=$dataSumName srcCopy
				CopyWaveAttributes(dataName,dimList[0],dataSumName,0)
				CopyWaveAttributes(dataName,dimList[1],dataSumName,1)
				break
			case 4:
				//Make/O/N=(DimSize(srcData,0)) temp1d
				DeletePoints dim,1, dimList

				if (!wholeRange)
					if (dim==0)
						duplicate /O /RMD=[start,finish][][][] srcData,srcCopy
					elseif (dim==1)
						duplicate /O /RMD=[][start,finish][][] srcData,srcCopy
					elseif (dim==2)
						duplicate /O /RMD=[][][start,finish][] srcData,srcCopy
					elseif (dim==3)
						duplicate /O /RMD=[][][][start,finish] srcData,srcCopy
					endif
				endif
				SumDimension /D=(dim)/DEST=$dataSumName srcCopy
				CopyWaveAttributes(dataName,dimList[0],dataSumName,0)
				CopyWaveAttributes(dataName,dimList[1],dataSumName,1)
				CopyWaveAttributes(dataName,dimList[2],dataSumName,2)
				break
		endswitch
	endif
End Function


//Copies wave attributes (wave units and wave scaling) from dimension OriginDim of OriginDataSet
//to the dimension DestinationDim of DestinationDataSet
Function CopyWaveAttributes(OriginDataSet,OriginDim,DestinationDataSet,DestinationDim)
	String OriginDataSet,DestinationDataSet
	Variable OriginDim,DestinationDim
	//first, check validity of the input
	if(WaveExists($OriginDataSet)&&WaveExists($DestinationDataSet))//check if input waves exists
		if((OriginDim<WaveDims($OriginDataSet))&&(DestinationDim<WaveDims($DestinationDataSet)))//and check if the given dimensions fit into the size of the wave
			//then, proceed copying the wave scaling and wave units
			Variable Offset=DimOffset($OriginDataSet,OriginDim)
			Variable Delta=DimDelta($OriginDataSet,OriginDim)
			String Units=WaveUnits($OriginDataSet,OriginDim)
			if(DestinationDim==0)
				SetScale/P x,Offset,Delta,Units,$DestinationDataSet
			elseif(DestinationDim==1)
				SetScale/P y,Offset,Delta,Units,$DestinationDataSet
			elseif(DestinationDim==2)
				SetScale/P z,Offset,Delta,Units,$DestinationDataSet
			elseif(DestinationDim==3)
				SetScale/P t,Offset,Delta,Units,$DestinationDataSet
			elseif(DestinationDim==-1)
				SetScale/I d,0,0,Units,$DestinationDataSet
			endif
		endif
	endif
End

////////////////////////////////////////////////////////////////////
// Normalize polar scan in a form of 2D or 3D wave E(em. angle)
////////////////////////////////////////////////////////////////////
// backNrPoints - number of top points used to calculate background
// backFactor -> data=data-background/(backFactor/100)
// normRange -> 100 means the whole energy range is used to calculate the integral
// normFactor -> data=data/integral^(normFactor/100)
function NormalizePolar(dataName,backNrPoints,backFactor,normRange,normFactor)
	String dataName
	Variable backNrPoints,backFactor,normFactor,normRange

	SetDataFolder ksDataFolder
	if(WaveExists($dataName))
		Variable i,j,k,background
		Variable	lastIndex,integralRange
		WAVE srcData = $dataName
		// Check if wave is 2 dimensional
		if(WaveDims($dataName)==2)
			Make/O/N=(DimSize(srcData,0)) data1d
			lastIndex=DimSize(srcData,0)-1
			integralRange = round(lastIndex/normRange)
			for (j=0;j<DimSize(srcData,1);j+=1)
				// background subtraction (top values - highest E kin)
				if (backFactor>0)
					data1d = srcData[p][j]
					srcData[][j] -= mean(data1d,lastIndex-backNrPoints,lastIndex)*(backFactor/100)
				endif
				// Normalize by integral
				if (normFactor>0)
					data1d = srcData[p][j]
					srcData[][j] /= mean(data1d,0,lastIndex-integralRange)^(normFactor/100)
				endif
			endfor
			KillWaves data1d
			// 3D map normalization
		elseif(WaveDims($dataName)==3)
			Make/O/N=(DimSize(srcData,2)) data1d
			lastIndex=DimSize(srcData,2)-1
			integralRange = round(lastIndex/normRange)
			// Work on source files
			String orgName = "SMPM"+dataName[3,strSearch(dataName,"T",0)-1]
			String list = WaveList(orgName+"*",";","")
			Wave/T listWv = ListToTextWave(list,";")    // Convert to text wave
			for (i=0;i<DimSize(listWv,0);i+=1)
				Wave data2D = $listWv[i]
				for (j=0;j<DimSize(data2D,1);j+=1)
					// background subtraction (top values - hjghest E kjn)
					if (backFactor>0)
						data1d[] = data2D[p][j]
						data2D[][j] -= mean(data1d,lastindex-backNrPoints,lastindex)*(backFactor/100)
					endif
					// Normalize by integral
					if (normFactor>0)
						data1d[]= data2D[p][j]
						data2D[][j] /= mean(data1d,0,lastindex-integralRange)^(normFactor/100)
					endif
				endfor
			endfor
			Make3Dmap(orgName)
			KillWaves data1d
		endif
	endif
end function

// Divide data by detector image
function NormalizeBgrImage(dataName,bgrName)
	String dataName,bgrName

	if (waveExists($dataName) && waveExists($bgrName))
		wave srcData = $dataName
		wave bgrData = $bgrName

		Duplicate/O bgrData, bgrDataRS
		Smooth 1, bgrDataRS
		Variable srcEDim = WaveDims(srcData)==2 ? 0 : 2
		if (DimSize(srcData,srcEDim)!=DimSize(bgrData,0) || DimSize(srcData,1)!=DimSize(bgrData,1))
			print "Dimension mismatch. \nIt means that the energy/angular channel number of the detector was changed."
			print "Data dimension: "+num2str(DimSize(srcData,srcEDim)) + " / " + num2str(DimSize(srcData,1))
			print "Background wave dimension: "+num2str(DimSize(bgrData,0)) + " / " + num2str(DimSize(bgrData,1))
			DoAlert 1, "Background data has different dimension size. \nDo you want to redimension it?"
			if (V_flag==2)	// User clicked NO
				abort
			endif
			Redimension/N=(DimSize(srcData,srcEDim),DimSize(srcData,1)) bgrDataRS
		endif
		if(WaveDims(srcData)==2)    // apply to single 2D scan
			srcData[][] /= bgrDataRS[p][q]
		elseif (WaveDims(srcData)==3)// apply to 3D data
			// Work on source files
			String orgName = "SMPM"+dataName[3,strSearch(dataName,"T",0)-1]
			String list = WaveList(orgName+"*",";","")
			Wave/T listWv = ListToTextWave(list,";")    // Convert to text wave
			Variable i
			for (i=0;i<DimSize(listWv,0);i+=1)
				Wave data2D = $listWv[i]
				data2D[][] /= bgrDataRS[p][q]
			endfor
			Make3Dmap(orgName)
		endif
	else
		print "Data not found"
	endif
	//KillWaves bgrDataRS
end function


// Divide data by detector image
function NormalizeBgr1D(dataName,bgrName)
	String dataName,bgrName

	if (waveExists($dataName) && waveExists($bgrName))
		wave srcData = $dataName
		wave bgrData = $bgrName

		ReduceDim(bgrName,0,"bgrData1D")
		wave bgrData1D = $"bgrData1D"
		bgrData1D /= DimSize(bgrData,0)

		Variable srcEDim = WaveDims(srcData)==2 ? 0 : 2
		if (DimSize(srcData,1)!=DimSize(bgrData1D,1))
			print "Dimension mismatch. The number of angular channels in bgr wave will be changed."
			Redimension/N=(DimSize(srcData,1)) bgrData1D
		endif
		if(WaveDims(srcData)==2)    // apply to single 2D scan
			srcData[][] /= bgrData1D[q]
		elseif (WaveDims(srcData)==3)// apply to 3D data
			// Work on source files
			String orgName = "SMPM"+dataName[3,strSearch(dataName,"T",0)-1]
			String list = WaveList(orgName+"*",";","")
			Wave/T listWv = ListToTextWave(list,";")    // Convert to text wave
			Variable i
			for (i=0;i<DimSize(listWv,0);i+=1)
				Wave data2D = $listWv[i]
				data2D[][] /= bgrData1D[q]
			endfor
			Make3Dmap(orgName)
		endif
	else
		print "Data not found"
	endif
	//KillWaves bgrDataRS
end function

/////////////////////////////////////
/// creates a series of 3D data sets
Function/S Make3Dmap(dataName)
	String dataName
	/////////////////////////
	SetDataFolder ksDataFolder
	String list=WaveList(dataName+"*",";","")		// Get the names of all waves
	/////////////////////////
	Make/O/N=(2,itemsInList(list)) goniometerPT     // prepare 2 rows wave for phi and theta values
	/////////////////////////
	Wave/T txtAttr=$ksTextAttributes
	/////////////////////////
	Variable i,j,k,l
	Variable numOfAreas=1                   // Numer of measured areas
	Variable numOfP=1                       // how many P points in one area
	Variable numOfE,numOfNu                 // number of Energy/detector angles points
	/////////////////////////
	// Determine the type of measurement
	string zone = "line"
	String stemp
	stemp=StringFromList(0,list)
	variable t1,t2
	t1 = str2num(txtAttr[%$stemp][%$ksaT])
	stemp=StringFromList(1,list)
	t2 = str2num(txtAttr[%$stemp][%$ksaT])
	if (abs(t1-t2)<0.01)
		zone = "region"
	endif
	print zone
	// Count the number of files and get the goniometer data
	for(i=0;i<itemsInList(list);i+=1)	// for each file (wave) do
		stemp=StringFromList(i,list)
		goniometerPT[0][i]=str2num(txtAttr[%$stemp][%$ksaP])
		goniometerPT[1][i]=str2num(txtAttr[%$stemp][%$ksaT])
		// check how many areas were measured
		if (cmpstr(zone,"region")==0)
			if (i>=1)
				if (abs(goniometerPT[1][i]-goniometerPT[1][i-1])>1)		// detect if theta steps are bigger that 1 degree
					numOfAreas+=1
				else
					if(numOfAreas==1)									// count the number of points in a single area
						numOfP+=1
					endif
				endif
			endif
		elseif (cmpstr(zone,"line")==0)
			numOfP=itemsInList(list)
		endif
	endfor
	//////////////////////////////////
	// Prepare 3D matrix with the data
	dataName=StringFromList(0,list)

	Variable deltaP, startP, startE, deltaE, deltaNu, startNu
	startP = goniometerPT[0][0]
	deltaP = goniometerPT[0][1]-startP
	deltaE = DimDelta($dataName,0)
	startE = DimOffset($dataName,0)
	numOfE = DimSize($dataName,0)
	numOfNu=DimSize($dataName,1)
	deltaNu=DimDelta($dataName,1)
	startNu=DimOffset($dataName,1)

	print "Number of: P point / areas / E / Nu"
	print numOfP, numOfAreas, numOfE, numOfNu
	if (cmpstr(zone,"region")==0)
		for( i=0; i<numOfAreas; i+=1)							// For each area make a Ang3D matrix
			Make/O/N=(numOfP,numOfNu,numOfE) Ang3D
			Ang3D[][][]=NaN
			SetScale/P x, startP, deltaP, "P (deg)" , Ang3D
			SetScale/P y, startNu,deltaNu,"nu (deg)", Ang3D
			SetScale/P z, startE, deltaE,"Ek (eV)",   Ang3D
			for (j=0;j<numOfP;j+=1)
				if (strlen(StringFromList(j+i*numOfP,list))>0)
					dataName=StringFromList(j+i*numOfP,list)
					Wave Im=$dataName
					// the symbols p, q, r and s take on the value of the row, column, layer and
					//chunk, respectively, of the element in the destination for which a value is being calculated
					Ang3D[j][][]=Im[r][q]
				endif
			endfor
			dataName=StringFromList(0, list)
			dataName="D3_"+dataName[4,strsearch(dataName,"_",0)]+"T"+num2str(Round(abs(goniometerPT[1][i*numOfP+1])))
			duplicate/O Ang3D $dataName
			print "### Region name: " + dataName
			if (i<numOfAreas-1)
				//DisplayECmap(dataName)
				print "To show this data use: SPEM->Browse loaded 3D maps \n"
			endif
			//DisplayECmap(dataName)
		endfor
	endif
	if (cmpstr(zone,"line")==0)
		Make/O/N=(itemsInList(list)) goniometerT
		goniometerT[]=goniometerPT[1][p]

		WaveStats goniometerT
		Variable Tav = V_avg
		Variable Tmin= V_min
		Variable Tmax= V_max
		// the size in y direction detector + number of T points
		Variable numOfNuT = numOfNu+ceil((Tmax-Tmin)/deltaNu)
		// move by the start Nu T value
		Variable startNuT =startNu+Tmin-Tav

		Make/O/N=(numOfP,numOfNuT,numOfE) Ang3D
		Ang3D[][][] = NaN
		SetScale/P x, startP,  deltaP, "P (deg)", Ang3D
		SetScale/P y, startNu+Tmin, deltaNu, "T + nu (deg)", Ang3D
		SetScale/P z, startE,  deltaE, Ang3D

		for (j=0;j<numOfP;j+=1)
			dataName=StringFromList(j,list)
			Wave Im=$dataName
			for(l=0;l<numOfNu; l+=1)
				Ang3D[j][l+ceil((goniometerT[j]-Tmin)/deltaNu)][] = Im[r][l]
			endfor
		endfor
		dataName=StringFromList(0, list)
		dataName="D3_"+dataName[4,strsearch(dataName,"_",0)]+"T"+num2str(Round(abs(Tav)))
		duplicate/O Ang3D $dataName
	endif
	KillWaves Ang3D
	return dataName
end
///////////////////////////////
