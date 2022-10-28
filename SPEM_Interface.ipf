#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3 // Use modern global access method and strict wave access.

// This file contains functions used for
// creation of dialog windows and handling the events

Menu "SPEM"
	help={"This menu gives you access to various SPEM procedures"}
	"Load data/F2",LoadHDFdata();
	"--------------------------------",;
	"Load sequence of data/F3",LoadHDFsequence();
	"--------------------------------",;
	"Set global value for Ef, offP and offT",ControlGLobalValues();
	"--------------------------------",;
	"Browse loaded Spectro Images", ControlBrowseSpectroImage();
	"Browse loaded 2D maps",ControlBrowse2Dmap();
	"Browse loaded 3D maps",ControlBrowse3Dmap();
	"--------------------------------",;
	"Show attributes",Edit :SPMData:txtAttributes.ld;
End

function DisplayByType(DatasetName,Type,WindowName)
	String DatasetName,Type,WindowName
	strswitch(Type)
		case ksSpectrum:
			Display1Dspectrum(dataSetName)
			break
		case ksImage:
			DisplaySpectroImage(DatasetName)
			break
		case ksImage4D:
			DisplaySpectroImage(DataSetName)
			break
		case ksPolar:
			Display2Dmap(DatasetName)
			break
		case ksAngularScan:
			//if (WaveDims($DataSetName)==3)
			String dataName = Make3DMap(dataSetName)
			DisplayECmap(dataName)
			//else
			//Display2Dmap(dataSetName)
			//endif
			break
	endswitch
end function

//////////////////////////////////////////////////
// 3D Make map control
// This function shows the dialog that allows to choose the data to create 3D map
Function ControlMake3Dmap()
	String nameSelected // Name of the selected data
	string nameListAll = WaveList("SMPM*",";","") // Prepare the list of the data
	String temp,nameList = ""
	// Remove items that form a sequence
	Variable i
	for(i=0;i<itemsInList(nameListAll);i+=1)
		temp = StringFromList(i,nameListAll,";")
		temp = StringFromList(0,temp,"_")
		if (WhichListItem(temp,nameList)==-1)
			nameList = addListItem(temp,nameList)
		endif
	endfor
	// Dialog window
	Prompt nameSelected "Select data", popup, nameList
	Doprompt "Select for transformation", nameSelected
	if (V_flag)							// User canceled
		abort
	endif
	Make3Dmap(nameSelected)			// Display the data
end

//////////////////////////////////////////////////
// 2D map Display dialog
// This function shows the dialog that allows to choose polar scan (2D) to display
Function ControlBrowse2Dmap()
	String nameSelected
	string nameListAll = WaveList("SMP*",";","") // Prepare the list of the data
	String temp,nameList = ""
	// Remove items that form a sequence
	Variable i
	for(i=0;i<itemsInList(nameListAll);i+=1)
		temp = StringFromList(i,nameListAll,";")
		if (StrSearch(temp,"SMPM",0)==-1)
			nameList = addListItem(temp,nameList)
		endif
	endfor
	// Dialog window
	Prompt nameSelected "Select 2D data", popup, nameList
	Doprompt "Select polar scan data to display", nameSelected
	if (V_flag)							// User canceled
		abort
	endif
	Display2Dmap(nameSelected)			// Display the data
end


//////////////////////////////////////////////////
// 3D map Display dialog
// The same as 2D map display but allows to select 3D data
Function ControlBrowse3Dmap()
	String nameSelected
	string nameListAll=WaveList("D3*",";","")		// Prepare the list of data
	String temp,nameList = ""
	// Remove items that form a sequence
	Variable i
	for(i=0;i<itemsInList(nameListAll);i+=1)
		temp = StringFromList(i,nameListAll,";")
		if (StrSearch(temp,"EK",0)==-1 && StrSearch(temp,"ld",0)==-1)
			nameList = addListItem(temp,nameList)
		endif
	endfor
	// Dialog window
	Prompt nameSelected "Select 3D data", popup, nameList
	Doprompt "Select 3D dataset to show",nameSelected
	if (V_flag)							// User canceled
		abort
	endif
	DisplayECmap(nameSelected)			// Display the data
end


///////////////////////////////////////////
// Display 1D spectrum
Function Display1Dspectrum(dataName)
	String dataName
	if (WaveExists($dataName))
		// Check the type of data in opened window
		String wName = StringFromList(0,TraceNameList("",";",1))
		if (strlen(wName)>5)
			WAVE/T atts=$ksTextAttributes
			String dataType = atts[%$wName][%$ksAcqType]
			if (stringMatch(dataType,ksSpectrum))
				DoAlert 1, "Do you want to append data to the current graph?"
				if (V_flag==2)	// User clicked NO
					Display /W=(200,50,450,250)
				endif
				AppendToGraph $dataName
			endif
		else
			Display /W=(200,50,450,250)
			AppendToGraph $dataName
			DelayUpdate
		endif
	endif
End function


///////////////////////////////////////////
// Display 2D ARPES data (polar scan)
Function Display2Dmap(dataName)
	String dataName
	String fldrSav=GetDataFolder(1)
	SetDataFolder ksDataFolder
	if(WaveExists($dataName))
		Display /W=(50,50,260,300)
		AppendImage $dataName;DelayUpdate
		SetAxis/A bottom;DelayUpdate
		SetAxis/A left
		ModifyGraph swapXY=1
		ModifyGraph mirror=1,standoff=0
		ModifyImage $dataName ctab= {*,*,$ksColorMap,0}
		addLegend(dataName) // add the legend to the plot
		ControlBar 30       // create controls
		Button btnFind,pos={10,5},size={40,20},title="Del",proc=eventDeleteData, help={"Delete the data from Igor"}
		Button btnNorm,pos={60,5},size={40,20},title="Norm",proc=eventNormPanel,  help={"Open normalization control panel"}
		Button btnShowInf,pos={110,5},size={40,20},title="Info",proc=eventShowInfo,  help={"Show info about the measurement"}
		Button btnTr2k,pos={160,5},size={40,20},title="To K",proc=eventTransform2K,  help={"Transform angles to k space"}
	endif
	SetDataFolder fldrSav
End function


///////////////////////////////////////////
// Display constant energy map Ek map vs emission angles
Function DisplayECmap(Wname) : Graph
	string Wname
	variable Ek,Ekin,Emin,Emax,KE,deltaE
	PauseUpdate; Silent 1
	Display /W=(50,50,405,380)
	AppendImage $WName
	string graphname=Winname(0,1)

	KE=DimSize($WName,2)
	Emin=DimOffset($WName,2)
	deltaE=DimDelta($WName,2)
	Emax=Emin+deltaE*(KE-1)
	ControlBar 50   // create controls
	Slider Ek,pos={10,5},size={135,50},fSize=8,proc=eventEKslider
	Slider Ek,limits={Emin,Emax,0},  value=22,vert=0,ticks=10

	ModifyGraph mirror=1,standoff=0
	ModifyImage $WName ctab= {*,*,$ksColorMap,0}
	//ModifyGraph width=280

	// Check if the data are in K space or vs emission angles
	if (StrSearch(WName,"_K",0)==-1)
		ValDisplay Ekin,pos={150,7},size={100,20},title="Ek (eV)"
		Button btnNorm,pos={150,25},size={50,20},title="Norm", proc=eventNormPanel, help={"Open normalization control panel"}
		Button btnShowInf,pos={205,25},size={50,20},title="Info",proc=eventShowInfo3D,  help={"Show info about the measurement"}
		Button btnBrowseK,pos={260,25},size={90,20},title="Browse E(ang)", proc=eventBrowseEK
		Button btnTrans,pos={260,5},size={90,20},title="Trans. to K", proc=eventTransform2K,help={"Transform angles to K space"}
		ValDisplay Ekin,pos={150,5},size={100,20},title="Ekin (eV)"
		Button btnShowInf,pos={200,25},size={55,20},title="Info",proc=eventShowInfo3D,  help={"Show info about the measurement"}
		ValDisplay Ekin,limits={500,500,500},barmisc={0,1000},value= _NUM:0.000
	else
		ValDisplay Ekin,pos={150,7},size={100,20},title="Ef (eV)"
		Button btnShowInf,pos={205,25},size={50,20},title="Info",proc=eventShowInfo3D,  help={"Show info about the measurement"}
		Button btnBrowseK,pos={260,25},size={90,20},title="Browse E(k)", proc=eventBrowseEK
		Button btnExtractImage,pos={260,5},size={90,20},title="Make plot",proc=eventPlotKxKy
		ValDisplay Ekin,pos={150,5},size={100,20},title="E - Ef (eV)"
		Button btnShowInf,pos={200,25},size={55,20},title="Info",proc=eventShowInfo3D,  help={"Show info about the measurement"}
		ValDisplay Ekin,limits={500,500,500},barmisc={0,1000},value= _NUM:0.000
	endif
	// make legend
	String orgName = "SMPM"+WName[3,strsearch(wName,"_",3)-1]+"_001"
	addLegend(orgName)
end


//////////////////////////////
// Windows events
//////////////////////////////
// Event transform to K button
function eventTransform2K(ctrlName): ButtonControl
	string ctrlName
	String srcName = StringFromList(0,ImageNameList("",";"))
	if (WaveDims($srcName)==3)
		TransformMerge2K(srcName)
	elseif (WaveDims($srcName)==2)
		Transform2Dscan2K(srcName)
	endif
	return 0
end

//////////////////////////////
// Event print the info button
function eventShowInfo(ctrlName): ButtonControl
	string ctrlName
	String srcName = StringFromList(0,ImageNameList("",";"))
	String dataName
	Variable ind
	ind = StrSearch(srcName,"sum",0)
	if (ind>0)
		dataName = srcName[0,ind-2]
	else
		dataName = srcName
	endif
	WAVE/T atts=$ksTextAttributes
	Print("<<< File info >>>")
	Print(atts[%$dataName][%$"Acquisition Type"])
	Print("P/T: "+atts[%$dataName][%$"aP"]+"° / "+atts[%$dataName][%$"aT"]+"°\r")
	Print("x/y/x: \t\t"+atts[%$dataName][%$"X"]+" / "+atts[%$dataName][%$"Y"]+" / "+atts[%$dataName][%$"Z"])
	Print("Ek start (ev): \t\t"+atts[%$dataName][%$"Offset 0"])
	Print("Pass Energy (ev): \t"+atts[%$dataName][%$"PE"])
	Print("Dwell time (s): \t"+atts[%$dataName][%$"Dwell"])
	Print("Scans: \t\t"+atts[%$dataName][%$"Scans"])
	Print("Start time: \t"+atts[%$dataName][%$"Date time start"])
	return 0
end function
//////////////////////////////
// Event print the info button for the 3D map
function eventShowInfo3D(ctrlName): ButtonControl
	string ctrlName
	String srcName = StringFromList(0,ImageNameList("",";"))
	String dataName
	Variable ind
	ind = StrSearch(srcName,"T",0)
	dataName = "SMPM" + srcName[3,ind-2] + "_001"
	String dataName2 = "SMPM" + srcName[3,ind-2] + "_002"
	WAVE/T atts=$ksTextAttributes
	Print("<<< File info >>>")
	Print("3D scan")
	Print("P/T start: "+atts[%$dataName][%$"aP"]+"° / "+atts[%$dataName][%$"aT"]+"°\r")
	Print("P delta : "+num2str(str2num(atts[%$dataName2][%$"aP"]) -str2num(atts[%$dataName][%$"aP"]))+"°\r")
	Print("x/y/x: \t\t"+atts[%$dataName][%$"X"]+" / "+atts[%$dataName][%$"Y"]+" / "+atts[%$dataName][%$"Z"])
	Print("Ek (ev): \t\t"+atts[%$dataName][%$"Offset 0"])
	Print("Pass Energy (ev): \t"+atts[%$dataName][%$"PE"])
	Print("Dwell time (s): \t"+atts[%$dataName][%$"Dwell"])
	Print("Scans: \t\t"+atts[%$dataName][%$"Scans"])
	Print("Start time: \t"+atts[%$dataName][%$"Date time start"])
	return 0
end function

// Event slider E kin moved
Function eventEKslider(slidername,E,event)
	String slidername
	Variable E, event
	String dataName
	Variable startE,deltaE,KE,pl,Ekind
	dataName=StringFromList(0,ImageNameList("",""))	// get data name
	KE=DimSize($dataName,2)
	startE=DimOffset($dataName,2)
	deltaE=DimDelta($dataName,2)
	If (E<=startE)
		pl=0
	elseif (E>=startE+deltaE*(KE-1))
		pl=KE-1
	else
		pl=round((E-startE)/deltaE)
	endif
	Ekind=(E-startE)/deltaE
	ValDisplay Ekin value=_NUM:E
	ModifyImage $dataName, plane=pl
End

// Event Browse EK button clicked -- create new sub window
Function eventBrowseEK(ctrlName): ButtonControl
	string ctrlName
	string dataName, NamePlane, planeDef
	Variable X1,Y1,X2,Y2
	dataName=StringFromList(0,ImageNameList("",""))
	////////////////////////////////////
	ShowInfo
	String CursInfA=CsrInfo(A)
	String CursInfB=CsrInfo(B)
	if((strlen(CursInfA)==0)||(strlen(CursInfB)==0))
		Cursor /I A, $dataName, DimOffset($dataName,0),DimOffset($dataName,1)
		Cursor /I B, $dataName, DimOffset($dataName,0),DimOffset($dataName,1)+(DimSize($dataName,1)-1)*DimDelta($dataName,1)
		CursInfA=CsrInfo(A)
		CursInfB=CsrInfo(B)
	endif
	X1=NumberByKey("POINT", CursInfA)
	Y1=NumberByKey("YPOINT",CursInfA)
	X2=NumberByKey("POINT", CursInfB)
	Y2=NumberByKey("YPOINT",CursInfB)
	//////////////////////////////////////
	NamePlane=dataName+"_EK"
	String browserName = "BEK_"+WinName(0,1)

	String graphName=WinName(0,1);
	SetWindow $graphName hook(myHook)=eventCursorMoved, hookEvents=7 //event when cursor were moved

	If(strlen(StringFromList(0,WinList(browserName,";","WIN:1")))==0)
		// create window if not exist
		DoWindow/K browserName
		Display/W=(450,50,700,350)/K=1/N=$browserName
		ControlBar 30
		print dataName
		if (StrSearch(dataName,"_K",0))
			Button ExtractPlane, pos={10,5}, size={100,20}, title="Make E(k) plot", proc=eventPlotEk
			CheckBox chBox2ndDiff, pos={120,10}, size={100,20}, title="2nd derivative", proc=event2ndDiff
			//SetWindow $browserName, hook(BrowserWindowClosed)=BrowserWindowClosed
			SetActiveSubwindow $browserName
			planeDef=dummitrans(dataName,X1,Y1,X2,Y2)
			wave M_ExtractedSurface
			duplicate/O M_ExtractedSurface $NamePlane
			AppendImage $NamePlane
			ModifyImage $NamePlane ctab= {*,*,$ksColorMap,0}
			ModifyGraph swapXY=1
			ModifyGraph margin(top)=30
			ModifyGraph mirror=1,standoff=0
		else
			//DoWindow/F $NewName
			//PlaneDef=dummitrans(Name,X1,Y1,X2,Y2)
			//duplicate/O M_ExtractedSurface $NamePlane
			//		TitleBox/Z Kplane, pos={200,3}, size={30,20}, title=PlaneDef
			//
			//		TitlBox=PlaneDef
		endif
	else
		DoWindow/F $browserName
	endif
end
/////////////////////////////
// cursor moved event - used on constant energy map and spectro image
function eventCursorMoved(H_Struct)
	STRUCT WMWinHookStruct &H_Struct
	if(H_Struct.eventCode==7)
		String dataName = StringFromList(0,ImageNameList("",";"))
		if (StrSearch(dataName,"_sum",0)!=-1)
			cursorUpdateSB(dataName,"B_"+WinName(0,1))
		else
			cursorUpdateEK(dataName,"BEK_"+WinName(0,1))
		endif
		DoUpdate
	endif
end

//////////////////////////////
// Update spectrum browser for Spectro Image
function cursorUpdateSB(dataName, browserName)
	string dataName, browserName
	//SetActiveSubwindow $browserName
	Variable X1,Y1,X2,Y2
	String CursInfo=CsrInfo(A)
	X1 = NumberByKey("POINT",CursInfo)
	Y1 = NumberByKey("YPOINT",CursInfo)
	CursInfo=CsrInfo(B)
	X2 = NumberByKey("POINT",CursInfo)
	Y2 = NumberByKey("YPOINT",CursInfo)

	Variable xum = DimOffset($dataName,0)+DimDelta($dataName,0)*X1
	Variable yum = DimOffset($dataName,1)+DimDelta($dataName,1)*Y1
	String/G gLastCursorPos = num2str(xum)+";"+num2str(yum)

	String originDataName = dataName[0,StrSearch(dataName,"_sum",0)-1]
	String spectrumAName = originDataName+"_SA"
	String spectrumBName = originDataName+"_SB"
	wave originData = $originDataName
	nvar integSize = root:SPMData:gProbeArea
	variable i,j
	if (numtype(X1)==0 && numtype(Y1)==0)
		wave spectrum = $spectrumAName
		// sum 1, 5 or 9 of data points depending on the integSize
		if (WaveDims(originData)==3)
			if (integSize == 1)
				spectrum[] = originData[p][X1][Y1]
			else
				X1 = X1==0 ? 1:X1
				Y1 = Y1==0 ? 1:Y1
				X1 = X1>=DimSize(originData,1)-1 ? DimSize(originData,1)-2:X1
				Y1 = Y1>=DimSize(originData,2)-1 ? DimSize(originData,2)-2:Y1
				spectrum[]=0
				for (j=-1;j<=1;j+=1)
					for (i=-1;i<=1;i+=1)
						if (integSize==9 || i!=j)
							spectrum[] += originData[p][X1+i][Y1+j]
						endif
					endfor
				endfor
			endif
			// Browse 4D image
		elseif (WaveDims(originData)==4)
			spectrum[][] = originData[p+1][q][X1][Y1]
		endif
	endif
	// Cursor B
	if (numtype(X2)==0 && numtype(Y2)==0 && WaveDims(originData)==3)
		wave spectrum = $spectrumBName
		// sum 1, 5 or 9 of data points depending on the integSize
		if (integSize == 1)
			spectrum[] = originData[p][X2][Y2]
		else
			X2 = X2==0 ? 1:X2
			Y2 = Y2==0 ? 1:Y2
			X2 = X2>=DimSize(originData,1)-1 ? DimSize(originData,1)-2:X2
			Y2 = Y2>=DimSize(originData,2)-1 ? DimSize(originData,2)-2:Y2
			spectrum[]=0
			for (j=-1;j<=1;j+=1)
				for (i=-1;i<=1;i+=1)
					if (integSize==9 || i!=j)
						spectrum[] += originData[p][X2+i][Y2+j]
					endif
				endfor
			endfor
		endif
	endif
end

//////////////////////////
// update EK map
function cursorUpdateEK(dataName, browserName)
	string dataName, browserName
	Variable X1,X2, Y1, Y2
	string CursInfo, nameplane, PlaneDef
	CursInfo=CsrInfo(A)
	X1=NumberByKey("POINT",CursInfo)
	Y1=NumberByKey("YPOINT",CursInfo)
	CursInfo=CsrInfo(B)
	X2=NumberByKey("POINT",CursInfo)
	Y2=NumberByKey("YPOINT",CursInfo)

	String/G gLastCursorPosK = num2str(X1)+"_"+num2str(Y1)+"_" + num2str(X2)+"_"+num2str(Y2)
	NamePlane=dataName+"_EK"
	if (numtype(X1)==0 && numtype(X2)==0 && numtype(Y1)==0 && numtype(Y2)==0)
		PlaneDef = dummitrans(dataName,X1,Y1,X2,Y2)
		wave M_ExtractedSurface
		duplicate/O M_ExtractedSurface $NamePlane

		// Calculate 2nd diff
		ControlInfo /W=$browserName chBox2ndDiff
		if (V_Value)
			Calc2ndDiff(NamePlane)
		endif
		X1 = DimOffset($dataName,0)+X1*DimDelta($dataName,0)
		Y1 = DimOffset($dataName,1)+Y1*DimDelta($dataName,1)
		X2 = DimOffset($dataName,0)+X2*DimDelta($dataName,0)
		Y2 = DimOffset($dataName,1)+Y2*DimDelta($dataName,1)

		DrawAction delete
		SetDrawEnv xcoord= bottom,ycoord= left
		SetDrawEnv linefgc= (65280,0,0);DelayUpdate
		SetDrawEnv dash= 3;DelayUpdate
		SetDrawEnv arrow= 1
		DrawLine X1,Y1,X2,Y2
		// Store the data in separate wave that could be used as a trace
		String kxkyLine = dataName + "_ld"
		Make/O N=(2), $kxkyLine = {Y1,Y2}
		Variable xdel
		xdel = X2-X1
		if (xdel==0)
			xdel = 0.0001
		endif
		SetScale/P x, X1,xdel,$kxkyLine
	endif
	//DoWindow /F $browserName
end


Function/S dummitrans(name, X1,Y1, X2,Y2)
	String name
	variable  X1,Y1, X2,Y2
	Variable SizeX, SizeY, SizeE
	Variable Xdelta, Ydelta, Edelta, XdeltaN
	Variable Xstart, Ystart, Estart
	variable X1p, Y1p, X2p, Y2p
	Variable N2, E1, E2
	String PlaneDefinition
	Wave tem=$name

	SizeX=DimSize(Tem,0)
	SizeY=DimSize(Tem,1)
	SizeE=DimSize(Tem,2)
	Xdelta=DimDelta(Tem,0)
	Ydelta=DimDelta(Tem,1)
	Edelta=DimDelta(Tem,2)
	Xstart=DimOffset(Tem,0)
	Ystart=DimOffset(Tem,1)
	Estart=DimOffset(Tem,2)

	N2=sqrt((X1-X2)*(X1-X2)+(Y1-Y2)*(Y1-Y2))+1
	//print N2
	X1p=Xstart+Xdelta*X1
	Y1p=Ystart+Ydelta*Y1
	X2p=Xstart+Xdelta*X2
	Y2p=Ystart+Ydelta*Y2
	E2=Estart
	E1=Estart+Edelta*(SizeE-1)


	ImageTransform/X={SizeE,N2,X1p,Y1p,E2,X1p,Y1p,E1,X2p,Y2p,E1} extractSurface tem

	XdeltaN=sqrt((X1p-X2p)*(X1p-X2p)+(Y1p-Y2p)*(Y1p-Y2p))/N2
	wave M_ExtractedSurface
	SetScale/P y, 0, XdeltaN,"kII",  M_ExtractedSurface
	SetScale/P x, Estart, Edelta, "eV", M_ExtractedSurface

	string K1,K2,K3,K4
	K1=num2str(X1p)
	K1=K1[0,strsearch(K1,".",0)+3]
	K2=num2str(Y1p)
	K2=K2[0,strsearch(K2,".",0)+3]
	K3=num2str(X2p)
	K3=K3[0,strsearch(K3,".",0)+3]
	K4=num2str(Y2p)
	K4=K4[0,strsearch(K4,".",0)+3]
	PlaneDefinition="Kxa="+K1+"Kya="+K2+"Kxb="+K3+"Kyb="+K4
	return PlaneDefinition
end


//////////////////////////////
// Extract E(k) - create new plot
Function eventPlotEk(ctrlName): ButtonControl
	string ctrlName
	String dataName = StringFromList(0,ImageNameList("",""))
	plotEK(dataName)
end
//////////////////////////////
// Show data as 2nd derivative
Function event2ndDiff(CB_Struct) : CheckBoxControl
	STRUCT WMCheckboxAction &CB_Struct
	String dataName = StringFromList(0,ImageNameList("",""))
	if(CB_Struct.checked)
		Calc2ndDiff(dataName)
	else
		//cursorUpdateEK(dataName,WinName(0,1))
	endif
	return 0
end

////////////////////////////////////////////
// Create Kx Ky Plot for given Kinetic energy
Function eventPlotKxKy(ctrlName): ButtonControl
	string ctrlName
	Variable E
	ControlInfo Ek; E=V_Value // Get the value of kinetic energy from the control
	String dataName = StringFromList(0,ImageNameList("",""))
	plotKxKy(dataName,E)
end function



////////////////////// SPEM Images interface ///////////////////
// Image Display dialog
Function ControlBrowseSpectroImage()
	String nameSelected // Name of the selected data
	string nameListAll=WaveList("SMI*",";","")	// Prepare the list of data
	nameListAll+=WaveList("SMM*",";","")	// Prepare the list of data
	String nameList = ""
	String temp
	Variable i
	for(i=0;i<itemsInList(nameListAll);i+=1)
		temp = StringFromList(i,nameListAll,";")
		if (StrSearch(temp,"sum",0)==-1 && StrSearch(temp,"SB",0)==-1 && StrSearch(temp,"SA",0)==-1 )
			nameList = addListItem(temp,nameList)
		endif
	endfor
	// Dialog window
	Prompt nameSelected "Select data", popup, nameList
	Doprompt "Select spectro image data to show", nameSelected
	if (V_flag)					// User canceled
		abort
	endif
	DisplaySpectroImage(nameSelected)			// Show data
end

///////////////////////////////////////////
// Display SPEM Images
Function DisplaySpectroImage(dataName) : Graph
	String dataName
	String dataSum

	SetDataFolder ksDataFolder
	if(WaveExists($dataName))
		dataSum=ImageSum(dataName)		//reduce image dimension (sum over the energies)
		if(strlen(dataSum))
			Display /W=(50,50,300,300)
			AppendImage $dataSum;DelayUpdate

			SetAxis/A bottom;DelayUpdate
			SetAxis/A left
			Label bottom "x (µm)";DelayUpdate
			ModifyGraph tickUnit(bottom)=1,notation(bottom)=1
			Label left "z (µm)";DelayUpdate
			ModifyGraph tickUnit(left)=1,notation(left)=1
			ModifyGraph mirror=1,standoff=0
			ModifyGraph width=0,height={Plan,1,left,bottom}
			// make legend
			addLegend(dataName)
			/////////////////////////////////
			ControlBar 30 // create controls
			Button btDel,pos={10,5},size={50,20},title="Del",proc=eventDeleteData,help={"Remove this dataset completely form Igor"}
			Button btAXZ,pos={70,5},size={50,20},title="Add XZ",proc=eventAddXZ,help={"Add the XZ position of selected ARPES measurements"}
			Button btnShowInf,pos={130,5},size={50,20},title="Info",proc=eventShowInfo,  help={"Show info about the measurement"}
			Button btnSpectra,pos={190,5},size={50,20},title="Spectra",proc=eventBrowseSpectra,  help={"Browse spectra for selected point"}
		endif
	endif
End function

// Event browse spectra
function eventBrowseSpectra(ctrlName): ButtonControl
	string ctrlName
	String dataName   = StringFromList(0,ImageNameList("",""))

	String originName = dataName[0,strsearch(dataName,"_sum",0)-1]
	wave originData   = $originName

	String spectrumAName = originName + "_SA"
	String spectrumBName = originName + "_SB"
	if (WaveDims(originData)==3)
		Make /O/N=(DimSize(originData,0)) $spectrumAName
		wave spectrumA = $spectrumAName
		spectrumA[] = originData[p][0][0]
		CopyWaveAttributes(originName,0,spectrumAName,0)
		Duplicate/O $spectrumAName $spectrumBName
		wave spectrumB = $spectrumBName
	elseif (WaveDims(originData)==4)
		Make /O/N=(DimSize(originData,0)-1,DimSize(originData,1)) $spectrumAName
		wave spectrumA = $spectrumAName
		spectrumA[][] = originData[p+1][q][0][0]
		CopyWaveAttributes(originName,0,spectrumAName,0)
		CopyWaveAttributes(originName,1,spectrumAName,1)
	endif

	ShowInfo
	String graphName=WinName(0,1);
	SetWindow $graphName hook(myHook)=eventCursorMoved, hookEvents=7 //event when cursor were moved

	String sbName = "BSI_"+graphName
	If(strlen(StringFromList(0,WinList(sbName,";","WIN:1")))==0)
		// create browser window if does not exist
		Variable/G root:SPMData:gProbeArea = 1
		DoWindow/K sbName
		Display/W=(700,50,920,300)/K=1/N=$sbName
		ControlBar 70
		SetActiveSubwindow $sbName
		Variable ekLimL,ekLimH,angLimL,angLimH
		if (WaveDims(spectrumA)==1)
			// 3D spectroimage
			Button bSetProbeArea, pos={10,5},size={95,20},title="Probe area", proc=eventSetProbeArea,help={"Set global value for the integration size"}
			CheckBox cRangeActive,pos={115,10},size={95,15},title="Select range",proc=eventRangeActive,  help={"Integrate spectroimage over selected energy/angle range"}
			SetVariable svEkstart,pos={10,30},size={95,20},title="Ek Start"
			SetVariable svEkstop, pos={115,30},size={95,20},title="Ek Stop"
			// Add both spectra and selected region to the graph
			AppendToGraph /C=(65535,0,0) spectrumA
			AppendToGraph /C=(0,40000,0) spectrumB
			// Setup the range selectors
			ekLimL = DimOffset(spectrumA,0)
			ekLimH = DimOffset(spectrumA,0)+DimSize(spectrumA,0)*DimDelta(spectrumA,0)
			ekLimL = round(ekLimL*10)*0.1
			ekLimH = round(ekLimH*10)*0.1
			SetVariable svEkstart, limits={ekLimL,ekLimH,0.2}, value = _NUM:ekLimL, proc = eventSetSiRange
			SetVariable svEkstop,  limits={ekLimL,ekLimH,0.2}, value = _NUM:ekLimH, proc = eventSetSiRange
			// Draw region
			DrawAction delete
			SetDrawEnv xcoord= bottom,ycoord= left,linefgc= (0,0,65535),fillpat= 0, linethick= 2.00;DelayUpdate
			DrawRect ekLimL,0,ekLimH,10
		else
			// 4D spectroimage
			Button makePlot,pos={10,5},size={95,20},title="Save as new plot",proc=eventMakePlotFrom4D,help={"Plot this data in a new window"}
			CheckBox cRangeActive,pos={115,10},size={95,15},title="Select range",proc=eventRangeActive,  help={"Integrate spectroimage over selected energy/angle range"}
			SetVariable svEkstart,pos={10,30},size={95,20},title="Ek Start"
			SetVariable svEkstop, pos={115,30},size={95,20},title="Ek Stop"
			SetVariable svAngstart,pos={10,50},size={95,20},title="Ang. Start"
			SetVariable svAngstop, pos={115,50},size={95,20},title="Ang. Stop"
			// Setup the range selectors
			ekLimL = DimOffset(spectrumA,0)
			ekLimH = DimOffset(spectrumA,0)+DimSize(spectrumA,0)*DimDelta(spectrumA,0)
			ekLimL = round(ekLimL*10)*0.1
			ekLimH = round(ekLimH*10)*0.1
			SetVariable svEkstart,	limits={ekLimL,ekLimH,0.2}, value = _NUM:ekLimL, proc = eventSetSiRange
			SetVariable svEkstop,  limits={ekLimL,ekLimH,0.2}, value = _NUM:ekLimH, proc = eventSetSiRange

			angLimL = DimOffset(spectrumA,1)
			angLimH = DimOffset(spectrumA,1)+DimSize(spectrumA,1)*DimDelta(spectrumA,1)
			angLimL = round(angLimL*10)*0.1
			angLimH = round(angLimH*10)*0.1
			SetVariable svAngstart,	limits={angLimL,angLimH,0.5}, value = _NUM:angLimL, proc = eventSetSiRange
			SetVariable svAngstop,  limits={angLimL,angLimH,0.5}, value = _NUM:angLimH, proc = eventSetSiRange

			// Add 2d spectrum and range selection plots to the graph
			AppendImage spectrumA
			ModifyImage '' ctab= {*,*,$ksColorMap,0}
			ModifyGraph swapXY=1
			// Draw region
			DrawAction delete
			SetDrawEnv xcoord= bottom,ycoord= left,linefgc= (65535,0,0),fillpat= 0, linethick= 2.00;DelayUpdate
			DrawRect angLimL,ekLimL,angLimH,ekLimH

		endif
		//TitleBox xPosLabel, pos={90,5},size={30,15},title="X=",frame=0
		ModifyGraph mirror=1,standoff=0
	else
		DoWindow/F $sbName
	endif
	return 0
end function

// event change the size of integration area in spectrum browser
Function eventSetSiRange(SV_Struct) : SetVariableControl
	STRUCT WMSetVariableAction &SV_Struct

	String bsiName=WinName(0,1)		//Get the graph name
	String dataName = StringFromList(0,TraceNameList("","",1)+ImageNameList("",""))
	String originName = dataName[0,strsearch(dataName,"_SA",0)-1]
	String spectrumAName = originName + "_SA"
	String spectrumBName = originName + "_SB"
	String dataSumName = originName + "_sum"

	Variable startEind,stopEind,startAind,stopAind
	Variable ek1,ek2,ang1,ang2
	if (WaveDims($originName)==3)
		wave spectrumA = $spectrumAName
		wave spectrumB = $spectrumBName
		// Get the values from two value controls
		ControlInfo /W=$bsiName svEkStart
		ek1 =  V_Value
		ControlInfo /W=$bsiName svEkStop
		ek2 =  V_Value
		// Draw region
		DrawAction delete
		SetDrawEnv xcoord= bottom,ycoord= left,linefgc= (0,0,65535),fillpat= 0, linethick= 2.00;DelayUpdate
		DrawRect ek1,0,ek2,max(WaveMax(spectrumA),WaveMax(spectrumB))

		// Find start and stop values
		startEInd = round((min(ek1,ek2)-DimOffset(spectrumA,0))/DimDelta(spectrumA,0))
		stopEInd = round((max(ek1,ek2)-DimOffset(spectrumA,0)-DimDelta(spectrumA,0))/DimDelta(spectrumA,0))
		// Finally refresh the original spectroimage
		ControlInfo /W=$bsiName cRangeActive
		// If Range Active check box selected
		if (V_Value)
			ReduceDim(originName,0,dataSumName,start=startEind,finish=stopEind)
		endif
	elseif (WaveDims($originName)==4)
		wave spectrumA = $spectrumAName
		// Get the values from value controls
		ControlInfo /W=$bsiName svEkStart
		ek1 =  V_Value
		ControlInfo /W=$bsiName svEkStop
		ek2 =  V_Value
		ControlInfo /W=$bsiName svAngStart
		ang1 =  V_Value
		ControlInfo /W=$bsiName svAngStop
		ang2 =  V_Value
		// Draw region
		DrawAction delete
		SetDrawEnv xcoord= bottom,ycoord= left,linefgc= (65535,0,0),fillpat= 0, linethick= 2.00;DelayUpdate
		DrawRect ang1,ek1,ang2,ek2
		// Find start and stop values
		startEind = round((min(ek1,ek2)-DimOffset(spectrumA,0))/DimDelta(spectrumA,0))
		stopEind = round((max(ek1,ek2)-DimOffset(spectrumA,0)-DimDelta(spectrumA,0))/DimDelta(spectrumA,0))
		startAind = round((min(ang1,ang2)-DimOffset(spectrumA,1))/DimDelta(spectrumA,1))
		stopAind = round((max(ang1,ang2)-DimOffset(spectrumA,1)-DimDelta(spectrumA,1))/DimDelta(spectrumA,1))
		//// Finally refresh the original spectroimage
		ControlInfo /W=$bsiName cRangeActive
		//// If Range Active check box selected
		if (V_Value)
			String dataSumE = originName+"_sumE"
			ReduceDim(originName,0,dataSumE,start=startEind,finish=stopEind)
			ReduceDim(dataSumE,0,dataSumName,start=startAind,finish=stopAind)
			killwaves $dataSumE
		endif
	endif
	return 0
End

// Re/Set the kinetic energy range active
function eventRangeActive(cb) : CheckBoxControl
	STRUCT WMCheckboxAction& cb

	if (cb.checked)
		STRUCT WMSetVariableAction st
		eventSetSiRange(st)
	else
		String bsiName=WinName(0,1)		//Get the graph name
		String dataName = StringFromList(0,TraceNameList("","",1)+ImageNameList("",""))
		String originName = dataName[0,strsearch(dataName,"_SA",0)-1]
		String dataSumName = originName + "_sum"

		if (WaveDims($originName)==3)
			ReduceDim(originName,0,dataSumName)
		elseif (WaveDims($originName)==4)
			string datasume = originname+"_sume"
			reducedim(originname,0,datasume)
			reducedim(datasume,0,datasumname)
			killwaves $datasume
		endif
	endif
	return 0
end function

// event change the size of integration area in spectrum browser
function eventSetProbeArea(ctrlName) : ButtonControl
	String ctrlName
	String integSize
	Prompt integSize "Integration size (px)", popup, "1;5;9"
	Doprompt "Set global value for integration size", integSize
	if (V_flag==1)
		abort
	endif
	Variable/G root:SPMData:gProbeArea = str2num(integSize)
	return 0
end function


// Make separate plot form 4D image
function eventMakePlotFrom4D(ctrlName) : ButtonControl
	String ctrlName
	String dataName = StringFromList(0,ImageNameList("",""))
	print dataName

	svar cursPos = gLastCursorPos
	String xpos = num2str(round(str2num(StringFromList(0,cursPos))));
	String ypos = num2str(round(abs(str2num(StringFromList(1,cursPos)))));
	//variable ypos =
	String dstDataName = dataName+"_"+xpos+"_"+ypos;
	print "New wave created: " + dstDataName
	Duplicate/O $dataName,$dstDataName
	Display /W=(50,50,250,300)
	AppendImage $dstDataName;DelayUpdate
	SetAxis/A bottom;DelayUpdate
	SetAxis/A left
	ModifyGraph swapXY=1
	ModifyGraph mirror=1,standoff=0
	ModifyImage '' ctab= {*,*,$ksColorMap,0}
	return 0
end function

// Event Delete data Button
function eventDeleteData(ctrlName): ButtonControl
	string ctrlName
	// Dialog window to confirm
	DoAlert 1, "Do you want to remove selected data from Igor?"
	if (V_flag==2)	// User clicked NO
		abort
	endif
	String name=StringFromList(0,ImageNameList("",""))	// get the window name
	RemoveFromAtts(name)			//Remove from table

	KillWindow $WinName(0,1)	//Close top window
	KillWaves /Z $name          //Remove data
	// Remove source data
	if (StrSearch(name,"sum",0))
		name = name[0,10]

		KillWaves /Z $name
	endif
	return 0
end function

// Event WORKing on IT
function eventAddXZ(ctrlName): ButtonControl
	string ctrlName
	Variable i,j,k,l
	// prepare the basename list
	// TODO add 2D data
	string stemp,basenamelst=""
	string listAll= WaveList("SMP*",";","")
	for(i=0;i<ItemsInList(listAll);i+=1)
		stemp = StringFromList(i,listAll,";")
		stemp = StringFromList(0,stemp,"_")
		if (StringMatch(stemp,"SMPM*")==0)
			basenamelst = AddListItem(StringFromList(i,listAll,";"),basenamelst)
		elseif (WhichListItem(stemp,basenamelst)==-1)
			basenamelst=AddListItem(stemp,basenamelst)
		endif
	endfor
	/////////////////////////
	String srcName
	Prompt srcName "Base Wave", popup, basenamelst
	Doprompt "Select 3D data to extract X Z postion", srcName
	if (V_flag==1)
		abort
	endif
	String list=WaveList(srcName+"*",";","")		// Get the of first wave in a list
	srcName = StringFromList(0,list)
	/////////////////////////
	String dataName=StringFromList(0,ImageNameList("",""))	// get the window name

	// Get data and store as wave
	WAVE/T atts=$ksTextAttributes
	variable srcX,srcZ
	srcX = str2num(atts[%$srcName][%$"X"])
	srcZ = str2num(atts[%$srcName][%$"Z"])
	String xzWaveName = dataName+"_xz_"+srcName
	//Make /O/N=(1,2) $xzWaveName={{srcX},{srcZ}}
	//AppendToGraph $xzWaveName[][1] vs $xzWaveName[][0]
	Make /O $xzWaveName={srcZ}
	SetScale/P x srcX, 1, $xzWaveName
	AppendToGraph $xzWaveName
	ModifyGraph marker=19,msize=2,mrkThick=1,rgb=(0,0,65535)
	ModifyGraph mode($xzWaveName)=3
	return 0
end function

/////////////////////////////////////////
// This window allows to set the global values for Fermi Level and offP and offT angles.
Function ControlGlobalValues()
	nvar gEf27 = root:SPMData:gEf27
	nvar gEf74 = root:SPMData:gEf74
	nvar gT0 = root:SPMData:gT0
	nvar gP90 = root:SPMData:gP90
	variable Ef27 = gEf27
	variable Ef74 = gEf74
	variable T0 = gT0
	variable P90 = gP90
	prompt Ef27, "Ef for 27 eV"
	prompt Ef74, "Ef for 74 eV"
	prompt P90, "P for normal emission"
	prompt T0, "T for normal emission"
	Doprompt    "Set Fermi Level",Ef27,Ef74,P90,T0
	gEf27 = Ef27
	gEf74 = Ef74
	gT0 = T0
	gP90 = P90
End

/////////////////////////////////////////
// This window allows to control the data normalization.
Function ControlNormPanel()
	String dataName = StringFromList(0,ImageNameList("",";"))

	NewDataFolder/O root:backNorm
	setDataFolder ksDataFolder
	Dfref dfr = root:backNorm
	if (WaveDims($dataName)==2)
		Duplicate/O $dataName  dfr:$dataName
	elseif (WaveDims($dataName)==3)
		// Get the names of all waves from this 3D measurement
		String orgName = "SMPM"+dataName[3,strSearch(dataName,"T",0)-1]
		String list = WaveList(orgName+"*",";","")
		Wave/T listWv = ListToTextWave(list,";")    // Convert to text wave
		Variable i
		for (i=0;i<DimSize(listWv,0);i+=1)
			Duplicate/O $listWv[i] dfr:$listWv[i]
		endfor
	endif
	NewDataFolder/O/S root:SPMData:NormControl
	// Create global variables used by the control panel.
	Variable temp = NumVarOrDefault(":backNrPoints",5)
	Variable/G backNrPoints = temp
	temp = NumVarOrDefault(":backFactor",100)
	Variable/G backFactor = temp
	temp = NumVarOrDefault(":normRange",100)
	Variable/G normRange = temp
	temp = NumVarOrDefault(":normFactor",100)
	Variable/G normFactor = temp
	String stemp = StrVarOrDefault(":bgrName","NONE")
	String/G bgrName = stemp
	String/G selDataName = dataName
	SetDataFolder ksDataFolder
	// Check if panel already exists
	If(strlen(WinList("PanelNorm",";","WIN:64"))==0)
		// Create window
		PauseUpdate; Silent 1 // building the window...
		NewPanel/W=(800,50,1040,300)/K=1 /N=PanelNorm as "Data normalization"

		TitleBox titleDataName,  pos={10,5},size={170,20},title="Selected data: "+dataName, frame=0
		TabControl tb, tabLabel(0)="Integration",pos={10,30}, size={220,190}
		TabControl tb, tabLabel(1)="Detector",value = 0,proc=eventTabProc
		// Tab integration
		TitleBox tab0_titleSubtract,  pos={20,30+40},size={200,20},title="Subtract background over the FL", frame=0
		SetVariable tab0_backNrSetVer,pos={20,50+40},size={200,20},title="Background num. of points"
		SetVariable tab0_backSetVer,  pos={20,70+40},size={200,20},title="Background correction (%)"

		TitleBox tab0_titleIntagral,      pos={20,100+40},size={200,20},title="Normalize by constant integral", frame=0
		SetVariable tab0_normRangeSetVer, pos={20,120+40},size={200,20},title="Norm energy range (%)"
		SetVariable tab0_normSetVer,      pos={20,140+40},size={200,20},title="Integral norm (%)"

		SetVariable tab0_backNrSetVer,    limits={2,80,1},   value=root:SPMData:NormControl:backNrPoints
		SetVariable tab0_backSetVer,      limits={0,100,5},  value=root:SPMData:NormControl:backFactor
		SetVariable tab0_normRangeSetVer, limits={10,100,5}, value=root:SPMData:NormControl:normRange
		SetVariable tab0_normSetVer,      limits={0,100,5},  value=root:SPMData:NormControl:normFactor
		// Tab detector:
		svar bgrName = root:SPMData:NormControl:bgrName
		TitleBox tab1_label1,pos={20,30+40},size={200,20},title="Background wave:",frame=0,disable=1
		TitleBox tab1_selectedBgr,pos={20,50+40},size={200,20},title=bgrName, frame=1,disable=1

		Button tab1_btnTest,pos={20,80+40},    size={200,20},title="Change background wave",proc=eventSelectBgrWave,disable=1
		PopUpMenu tab1_normType,pos={20,120+40},size={200,20},title="Norm. type: ",value="by 2D image (SNAP);by 1D profile (SWEP scan)",disable=1

		Button btnDoIt,pos={10,220},    size={60,25},proc=eventNormCompute,title="Do it"
		Button btnSave,pos={90,220},    size={60,25},proc=eventNormSave,title="Save"
		Button btnRevert, pos={170,220},size={60,25},proc=eventRevertNorm,title="Revert"
	else
		SetActiveSubWindow panelNorm
	endif
End
Function eventSelectBgrWave(ctrlName): ButtonControl
	String ctrlName
	// Prepare list of SMP waves
	string nameListAll = WaveList("SMP*",";","") // Prepare the list of the data
	String nameSelected,temp,nameList = ""
	// Remove items that form a sequence
	Variable i
	for(i=0;i<itemsInList(nameListAll);i+=1)
		temp = StringFromList(i,nameListAll,";")
		if (StrSearch(temp,"SMPM",0)==-1)
			nameList = addListItem(temp,nameList)
		endif
	endfor
	// Dialog window
	Prompt nameSelected "Select wave for background correction", popup, nameList
	Doprompt "Select polar scan data to display", nameSelected
	if (V_flag)							// User canceled
		abort
	endif
	if (WaveDims($nameSelected)==2)
		svar bgrName = root:SPMData:NormControl:bgrName
		bgrName = nameSelected
		ModifyControl tab1_selectedBgr title=bgrName
	else
		print "Wrong data dimension! Background wave should be 2D."
	endif
end Function

Function eventTabProc(tca) : TabControl
	STRUCT WMTabControlAction &tca

	switch (tca.eventCode)
		case 2:                                  // Mouse up
			Variable tabNum = tca.tab           // Active tab number
			Variable isTab0 = tabNum==0
			Variable isTab1 = tabNum==1

			ModifyControl tab0_titleSubtract disable=!isTab0
			ModifyControl tab0_backNrSetVer disable=!isTab0
			ModifyControl tab0_backSetVer disable=!isTab0
			ModifyControl tab0_titleIntagral disable=!isTab0
			ModifyControl tab0_normRangeSetVer disable=!isTab0
			ModifyControl tab0_normSetVer disable=!isTab0
			ModifyControl tab0_backNrSetVer disable=!isTab0
			ModifyControl tab0_backSetVer disable=!isTab0
			ModifyControl tab0_normRangeSetVer disable=!isTab0
			ModifyControl tab0_normSetVer disable=!isTab0

			ModifyControl tab1_btnTest disable=!isTab1
			ModifyControl tab1_label1 disable=!isTab1
			ModifyControl tab1_selectedBgr disable=!isTab1
			ModifyControl tab1_normType disable=!isTab1
			break
	endswitch

	return 0
End

// Event norm panel
function eventNormPanel(ctrlName): ButtonControl
	string ctrlName
	// Dialog window to confirm
	ControlNormPanel()
	return 0
end function

// Compute normalized data for top window
Function eventNormCompute(ctrlName): ButtonControl
	String ctrlName
	svar dataName = root:SPMData:NormControl:selDataName
	// restore data from back copy
	Dfref dfr = root:backNorm
	if (WaveDims($dataName)==2)
		Duplicate/O dfr:$dataName  $dataName
	elseif (WaveDims($dataName)==3)
		String orgName = "SMPM"+dataName[3,strSearch(dataName,"T",0)-1]
		String list = WaveList(orgName+"*",";","")
		Wave/T listWv = ListToTextWave(list,";")    // Convert to text wave
		Variable i
		for (i=0;i<DimSize(listWv,0);i+=1)
			Duplicate/O dfr:$listWv[i] $listWv[i]
		endfor
	endif

	dfr = root:SPMData:NormControl
	nvar backNrPoints = dfr:backNrPoints
	nvar backFactor = dfr:backFactor
	nvar normRange = dfr:normRange
	nvar normFactor = dfr:normFactor
	svar bgrName = dfr:bgrName

	ControlInfo /W=panelNorm tb
	Variable tabNr = V_Value
	ControlInfo /W=panelNorm tab1_normType
	Variable normType = V_Value
	if (tabNr==0)
		print "Integral norm."
		NormalizePolar(dataName,backNrPoints,backFactor,normRange,normFactor)
	elseif(tabNr==1)
		print "Background correction"
		if (normType==1)
			NormalizeBgrImage(dataName,bgrName)
		elseif (normType==2)
			NormalizeBgr1D(dataName,bgrName)
		endif
	endif
End

// Compute normalized data for the whole set
Function eventNormSave(ctrlName): ButtonControl
	String ctrlName
	KillWindow panelNorm
	svar dataName = root:SPMData:NormControl:selDataName
	Dfref dfr = root:backNorm
	if (WaveDims($dataName)==2)
		Killwaves/Z dfr:$dataName
	elseif (WaveDims($dataName)==3)
		String orgName = "SMPM"+dataName[3,strSearch(dataName,"T",0)-1]
		String list = WaveList(orgName+"*",";","")
		Wave/T listWv = ListToTextWave(list,";")    // Convert to text wave
		Variable i
		for (i=0;i<DimSize(listWv,0);i+=1)
			KillWaves/Z dfr:$listWv[i]
		endfor
	endif
End

// Revert changes
Function eventRevertNorm(ctrlName): ButtonControl
	String ctrlName
	KillWindow panelNorm

	svar dataName = root:SPMData:NormControl:selDataName
	Dfref dfr = root:backNorm

	if (WaveDims($dataName)==2)
		Duplicate/O dfr:$dataName $dataName
		Killwaves/Z dfr:$dataName
	elseif (WaveDims($dataName)==3)
		String orgName = "SMPM"+dataName[3,strSearch(dataName,"T",0)-1]
		String list = WaveList(orgName+"*",";","")
		Wave/T listWv = ListToTextWave(list,";")    // Convert to text wave
		Variable i
		for (i=0;i<DimSize(listWv,0);i+=1)
			Duplicate/O dfr:$listWv[i] $listWv[i]
		endfor
		Make3Dmap(orgName)
	endif
End


// add legend to the graph
// the legend contains different information depending on the acquisition type
Function addLegend(dataName)
	String dataName
	String dataType

	String slegend
	WAVE/T atts=$ksTextAttributes

	dataType = atts[%$dataName][%$ksAcqType]

	strswitch(dataType)
		case ksSpectrum:
			break
		case ksPolar:
			slegend = "\Z08P/T = "
			slegend += num2str(round(str2num(atts[%$dataName][%$ksaP])*10)*0.1)+"° / "
			slegend += num2str(round(str2num(atts[%$dataName][%$ksaT])*10)*0.1)+"°\rx = "
			slegend += num2str(round(str2num(atts[%$dataName][%$ksX])*10)*0.1)+"\rz = "
			slegend += num2str(round(str2num(atts[%$dataName][%$ksZ])*10)*0.1)
			break
		case ksAngularScan:
			slegend = "\Z08x = "
			slegend += num2str(round(str2num(atts[%$dataName][%$ksX])*10)*0.1)+"\rz = "
			slegend += num2str(round(str2num(atts[%$dataName][%$ksZ])*10)*0.1)
			break
		case ksImage:
			slegend = "\Z08P/T = "
			slegend += num2str(round(str2num(atts[%$dataName][%$ksaP])*10)*0.1)+"° / "
			slegend += num2str(round(str2num(atts[%$dataName][%$ksaT])*10)*0.1)+"° \rEk = "
			slegend += num2str(round((	 str2num(atts[%$dataName][%$ksOffset0])  \
			+str2num(atts[%$dataName][%$ksDelta0])   \
			*dimSize($dataName,0)*0.5)*10)*0.1)+" eV"
			break
		case ksImage4D:
			slegend = "\Z08P/T = "
			slegend += num2str(round(str2num(atts[%$dataName][%$ksaP])*10)*0.1)+"° / "
			slegend += num2str(round(str2num(atts[%$dataName][%$ksaT])*10)*0.1)+"°"
			break
		default:
			slegend = ""
		endswitch
		TextBox/C/N=text0 /A=RB slegend

	end function


