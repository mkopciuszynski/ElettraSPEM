#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3   // Use modern global access method and strict wave access.

function plotEk(wName)
	String wName

	wave srcData = $wName
	String dstName = wName + "_plotEK"
	if (strsearch(wName,"D3",0)>=0)
		svar cursPos = gLastCursorPosK
		dstName += "_"+CursPos
	endif
	Duplicate/O srcData, $dstName

	wave data = $dstName

	MatrixOp/O data=ReplaceNaNs(data,0)
	// Setup the new plot EK
	PauseUpdate; Silent 1
	Display /W=(50,50,250,300)
	AppendImage data
	
	ModifyGraph tickUnit=1
	ModifyGraph swapXY=1
	Label left "\F'Arial' \Z10 E - E\BF \M (\Z09eV\Z10)"
	Label bottom "\F'Arial' \Z10 k\B∥ \M (\Z091/Å\M\Z10)"
	ModifyGraph mirror=1,standoff=0
	ModifyImage $dstName ctab= {*,*,$ksColorMap,0}
	//////////////////////////////////////////////////////////////////////////

end function


function plotKxKy(wName,E)
	String wName
	Variable E

	wave srcData = $wName
	Variable startE, deltaE, sizeE
	sizeE=DimSize(srcData,2); deltaE=DimDelta(srcData,2); startE=DimOffset(srcData,2)
	variable plane = round((E-startE)/deltaE)

	Make/O/N=(DimSize(srcData,0),DimSize(srcData,1)) dstMap
	// transfer one plane of 3D wave into 2D map
	dstMap[][] = srcData[p][q][plane]
	// set scale and offset of a new wave
	SetScale/P x, DimOffset(srcData,0), DimDelta(srcData,0), dstMap
	SetScale/P y, DimOffset(srcData,1), DimDelta(srcData,1), dstMap
	MatrixOp/O dstMap=ReplaceNaNs(dstMap,0)

	String dstName = wName +"_plotKK_"+ num2str(round(-E*1000))
	Duplicate/O dstMap, $dstName
	KillWaves dstMap

	//////////////////////////////////////////////////////////////////////////
	// Setup the new plot KXKY
	PauseUpdate; Silent 1 // building window...
	Display /W=(50,50,300,300)
	AppendImage $dstName
	Label left "\F'Arial' \Z10 k\By \M (\Z091/Å\M\Z10)"
	Label bottom "\F'Arial' \Z10 k\Bx \M (\Z091/Å\M\Z10)"
	ModifyGraph mirror=1,standoff=0
	ModifyImage $dstName ctab= {*,*,$ksColorMap,0}
	//ModifyGraph width={Aspect,1}
	ModifyGraph height={Plan,1,left,bottom}

	//////////////////////////////////////////////////////////////////////////
	// Add the line for extracted EK
	String kxkyLine = wName + "_ld"
	if (waveExists($kxkyLine))
		AppendToGraph $kxkyLine
	endif

	TextBox/C/N=textE/A=LB num2str(E)+" eV"

end function
