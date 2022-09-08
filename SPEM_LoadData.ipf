#pragma rtGlobals=3 // Use modern global access method.

/////////////////////////////////////////////////////////////
//This function, calls LoadHDF5Table() for a general load of data into Igor structure and displays it
function LoadHDFdata()
	//Create wave with text atributes if not exist
	if(!WaveExists($ksTextAttributes))
		InitAttTable()
	endif
	// Dialog to select the file to load
	variable container
	Open/D/R /F="HDF5 data file:.hdf5;" Container
	String filePathName = S_fileName
	// Open HDF5 file
	if (strlen(filePathName))
		String dataset = LoadHDF5Table(filePathName)
		if(strlen(dataset))
			WAVE/T Atts=$ksTextAttributes
			String Type=Atts[%$Dataset][%$ksAcqType]
			if (!stringMatch(type,ksAngularScan))
				DisplayByType(Dataset,Type,"")
			endif
		else
			print "The data was not loaded properly"
		endif
	endif
end

/////////////////////////////////////////////////////////////
// Load a sequence of files (3D angular scan)
function LoadHDFsequence()
	//Create wave with text atributes if not exist
	if(!WaveExists($ksTextAttributes))
		InitAttTable()
	endif
	String pathName,fileName,filePathName
	String fileList,dataSetName
	//////////////////////////////////////////
	Variable Container
	Open/D/R /F="HDF5 data file:.hdf5;" Container
	filePathName = S_fileName
	//////////////////////////////////////////
	pathName = ParseFilePath(1,filePathName,":",1,0)
	fileName = ParseFilePath(3,filePathName,":",0,0)
	dataSetName=fileName[0,strsearch(fileName,"_",0)-1]
	NewPath/O myPath, PathName
	fileList = IndexedFile(myPath,-1,".hdf5")
	//////////////////////////////////////////
	Variable i=0,numberOfFiles=0
	for (i=0;i<itemsInList(fileList);i+=1)
		fileName=StringFromList(i,Filelist)
		if(stringmatch(fileName,"SMPM*")!=1)
			FileList= RemoveFromList(fileName, FileList)
		endif
	endfor
	numberOfFiles = itemsInList(fileList)
	//////////////////////////////////////////
	if (numberOfFiles>=1)
		for(i=0;i<numberOfFiles;i+=1)
			fileName=StringFromList(i,fileList)
			filePathName = pathName+fileName
			LoadHDF5Table(filePathName)
		endfor
		DisplayByType(dataSetName,ksAngularScan,"")
	endif
	//////////////////////////////////////////
end Function

/////////////////////////////////////
// Open single HDF5 file
function/T LoadHDF5Table(filePathName)
	String filePathName
	Variable fileID,i,j,CurrentAtt=1,CurrentDts=-1
	String DataSetName=""
	String fldrSav=GetDataFolder(1)
	SetDataFolder ksDataFolder
	SetWaveLock 0, $ksTextAttributes //unlock textAttributes

	HDF5OpenFile/R/Z fileID as filePathName

	if (V_Flag==0)
		DataSetName=ParseFilePath(3,filePathName,":",0,0)
		WAVE/T txtAtt=$ksTextAttributes
		//check if the dataset with the given name already exists
		for(i=0;i<DimSize($ksTextAttributes,0);i+=1)
			if(cmpstr(txtAtt[i][0],DataSetName)==0)
				CurrentDts = i
			endif
		endfor
		Struct HDF5DataInfo di //strucutre to read into the dataset properties
		InitHDF5DataInfo(di)
		HDF5LoadData/O/Q/Z fileID, DataSetName //Read dataset from file
		if(V_flag==0)
			HDF5ListAttributes/Z fileID, DataSetName //Get the lists of all attributes of the loaded dataset
			if(V_Flag==0)
				if(CurrentDts==-1)
					CurrentDts=DimSize($ksTextAttributes,0)
					Redimension/N=((CurrentDts+1),-1) $ksTextAttributes
				endif
				txtAtt[CurrentDts][0]=DataSetName
				//For each attribute loads it's value(s)
				int HA
				String CurrentAttribute=""
				for(i=0;i<ItemsInList(S_HDF5ListAttributes);i+=1)
					CurrentAttribute=StringFromList(i,S_HDF5ListAttributes)
					HA=HDF5AttributeInfo(fileID,"/"+DataSetName,2,CurrentAttribute,0,di)	//Obtains type of the current attribute
					HDF5LoadData/A=CurrentAttribute/O/Q fileID, DataSetName             //Loads attribute value(s) into wave CurrentAttribute
					CurrentAtt=Attribute2tableNumber(CurrentAttribute)
					CurrentAttribute=CurrentAttribute[0,30]
					if(V_flag==0)
						// Check if attribute is string or number, and assign it's value with appropriate conversion
						strswitch(di.datatype_str)
							case "string":
								WAVE/T tmpstr=$CurrentAttribute//if string - creates text wave
								for(j=0;j<DimSize($CurrentAttribute,0);j+=1)
									txtAtt[CurrentDts][CurrentAtt]=tmpstr[j]
									CurrentAtt+=1
								endfor
								break
							default:
								WAVE tmp=$CurrentAttribute//if number - creates number wave
								for(j=0;j<DimSize($CurrentAttribute,0);j+=1)
									txtAtt[CurrentDts][CurrentAtt]=num2str(tmp[j])
									CurrentAtt+=1
								endfor
								break
						endswitch
						KillWaves/Z $CurrentAttribute
					endif
				endfor
			endif //HDFListAttributes
		endif //HDFLoadData
		WAVE tmp=$""
		WAVE/T tmpstr=$""
		//Sets a dataset name as a dimension label for the new row in the table
		SetDimLabel 0,CurrentDts,$DataSetName,$ksTextAttributes
		//Scale the axis
		SetScale/P x,str2num(txtAtt[CurrentDts][%$ksOffset0]),str2num(txtAtt[CurrentDts][%$ksDelta0]),txtAtt[CurrentDts][%$ksUnits0],$DatasetName
		if(strlen(txtAtt[CurrentDts][%$ksName1])!=0)
			SetScale/P y,str2num(txtAtt[CurrentDts][%$ksOffset1]),str2num(txtAtt[CurrentDts][%$ksDelta1]),txtAtt[CurrentDts][%$ksUnits1],$DatasetName
		endif
		if(strlen(txtAtt[CurrentDts][%$ksName2])!=0)
			SetScale/P z,str2num(txtAtt[CurrentDts][%$ksOffset2]),str2num(txtAtt[CurrentDts][%$ksDelta2]),txtAtt[CurrentDts][%$ksUnits2],$DatasetName
		endif
		if(strlen(txtAtt[CurrentDts][%$ksName3])!=0)
			SetScale/P t,str2num(txtAtt[CurrentDts][%$ksOffset3]),str2num(txtAtt[CurrentDts][%$ksDelta3]),txtAtt[CurrentDts][%$ksUnits3],$DatasetName
		endif
		if(cmpstr(txtAtt[CurrentDts][%$ksAcqType],ksImage)==0)
			DeletePoints/M=0 0,1,$DataSetName
		endif
		//if(cmpstr(txtAtt[CurrentDts][%$ksAcqType],ksImage4d)==0)   // 4D image treatment makes sum of angle and energy channels for each point and creates a new 2D wave with the image name_sum
		////Sum4D(DatasetName,str2num(txtAtt[CurrentDts][%$ksDwell]))
		//endif
		//Reduce extra-dimension of spectrum
		if(cmpstr( txtAtt[CurrentDts][%$ksAcqType],ksSpectrum)==0)
			Redimension/N=((DimSize($DatasetName,0)),0)  $DatasetName
		endif
		//Normalize analyzer counts per dwell and number of sweeps
		if(cmpstr(txtAtt[CurrentDts][%$ksDetector],"Analyser")==0)
			WAVE Dts=$DataSetName
			Dts/=str2num(txtAtt[CurrentDts][%$ksDwell])
			Dts/=str2num(txtAtt[CurrentDts][%$ksScans])
			WAVE Dts=$""
		endif
		//report success
		print  txtAtt[CurrentDts][%$ksAcqType],DataSetName, "loaded"
		WAVE/T txtAtt=$""
	elseif(V_flag!=-1)
		//The data load failed with V_flag error code
	endif
	SetWaveLock 1, $ksTextAttributes
	SetDataFolder fldrSav
	return DataSetName
end function

//
function Attribute2tableNumber(CurrentAt)
	string CurrentAt
	strswitch(CurrentAt)
		case "Sample ID":
			return 1
			break
		case "Acquisition Type":
			return 2
			break
		case "Dim0 Name Units":
			return 3
			break
		case "Dim0 Values":
			return 5
			break
		case "Dim1 Name Units":
			return 7
			break
		case "Dim1 Values":
			return 9
			break
		case "Dim2 Name Units":
			return 11
			break
		case "Dim2 Values":
			return 13
			break
		case "Stage Coord (XYZR)":
			return 15
			break
		case "Angular Coord":
			return 19
			break
		case "Detector":
			return 21
			break
		case "MCP Voltage":
			return 22
			break
		case "Lens Mode":
			return 23
			break
		case "Ep (eV)":
			return 24
			break
		case "Dwell Time (s)":
			return 25
			break
		case "N of Scans":
			return 26
			break
		case "DET Limits":
			return 27
			break
		case "Energy Window (eV)":
			return 29
			break
		case "Temperature (K)":
			return 31
			break
		case "Pressure (mbar)":
			return 33
			break
		case "Ring En (GeV) GAP (mm) Photon (eV)":
			return 34
			break
		case "Ring Current (mA)":
			return 37
			break
		case "Date Time Start Stop":
			return 39
			break
		case "Calibration":
			return 41
			break
		case "Defocus":
			return 42
			break
		case "Dim3 Name Units":
			return 43
			break
		case "Dim3 Values":
			return 45
			break
	endswitch
end

//Returns the value of Attribute corresponding to DatasetName
//Returns empty string, if now attribute with such a name is found
//Retuns the name of the dataset name if empty string is passed as an attribute
function/T GetByName(DatasetName,Attribute)
	String DatasetName,Attribute
	WAVE/T Atts=$ksTextAttributes
	Variable i,AttNumber=-1
	String Found=""
	String fldrSav=GetDataFolder(1)
	SetDataFolder ksSPMDataF
	if(strlen(Attribute))
		AttNumber=WhichListItem(Attribute,ksAtts)
		if(AttNumber==-1)
			WAVE/T Atts=$""
			SetDataFolder fldrSav
			return Found
		endif
	else
		AttNumber=0
	endif
	for(i=0;i<DimSize($ksTextAttributes,0);i+=1)
		if(cmpstr(Atts[i][0],DatasetName)==0)
			Found=Atts[i][AttNumber]
		endif
	endfor
	WAVE/T Atts=$""
	SetDataFolder fldrSav
	return Found
end


// NOR working properly!!!
//Deletes absent datasets from table
function VerifyTable()
	String fldrSav=GetDataFolder(1)
	SetDataFolder ksDataFolder
	String SMWaves=WaveList("SM*",";","DP:1")
	Variable i
	WAVE/T txtAtts=$ksTextAttributes
	SetWaveLock 0, $ksTextAttributes
	for(i=0;i<DimSize($ksTextAttributes,0);i+=1)
		if(WhichListItem(txtAtts[i][0],SMWaves)==-1)
			DeletePoints/M=0 i,1,$ksTextAttributes
			i-=1
		endif
	endfor
	WAVE/T txtAtts=$""
	SetWaveLock 1, $ksTextAttributes
	SetDataFolder fldrSav
end function

//Create wave with text attributes of each data file
function InitAttTable()
	String fldrSav=GetDataFolder(1)
	SetDataFolder ksSPMDataF
	Make/O/T/N=(0,ItemsInList(ksAtts)) $ksTextAttributes
	variable i
	for(i=0;i<ItemsInList(ksAtts);i+=1)
		SetDimLabel 1,i,$StringFromList(i,ksAtts),$ksTextAttributes
	endfor
	SetWaveLock 1, $ksTextAttributes
	SetDataFolder fldrSav
end function

//Kills old attribute table, if any, and creates the new one
function ResetAttTable()
	String fldrSav=GetDataFolder(1)
	SetDataFolder ksSPMDataF
	if(WaveExists($ksTextAttributes))
		SetWaveLock 0, $ksTextAttributes
	endif
	Make/O/T/N=(0,ItemsInList(ksAtts)) $ksTextAttributes
	variable i
	for(i=0;i<ItemsInList(ksAtts);i+=1)
		SetDimLabel 1,i,$StringFromList(i,ksAtts),$ksTextAttributes
	endfor
	SetWaveLock 1, $ksTextAttributes
	SetDataFolder fldrSav
end function

// Removes dataset from table
Function RemoveFromAtts(DataSetName)
	String DataSetName
	SetWaveLock 0, $ksTextAttributes
	Variable vReturn=0
	Variable ind
	ind = strsearch(DataSetName,"_",8)
	if (ind!=-1)
		DataSetName = DataSetName[0,ind-1]
	endif
	ind=FindDimLabel($ksTextAttributes,0,DataSetName)
	if(ind>=0)
		DeletePoints/M=0 ind,1,$ksTextAttributes
		vReturn = 1
	endif
	if (DimSize($ksTextAttributes,0)==0)
		KillWaves/Z $ksTextAttributes
	else
		SetWaveLock 1, $ksTextAttributes
	endif
End Function


//converts date and time stamp "DD/MM/YY hh:mm:ss" to the number of seconds since some date
function DateTime2secs(DateTimeStamp)
	String DateTimeStamp
	String myDate, myTime
	myDate=StringFromList(0,DateTimeStamp," ")
	myTime=StringFromList(1,DateTimeStamp," ")
	Variable hours,minutes,seconds
	hours=str2num(StringFromList(0,myTime,":"))
	minutes=str2num(StringFromList(1,myTime,":"))
	seconds=str2num(StringFromList(2,myTime,":"))
	Variable Day, Month,Year
	Day=str2num(StringFromList(0,myDate,"/"))
	Month=str2num(StringFromList(1,myDate,"/"))
	Year=str2num(StringFromList(2,myDate,"/"))
	return date2secs(Year,Month,Day)+3600*hours+60*minutes+seconds
End
