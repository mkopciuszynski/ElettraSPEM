#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3 // Use modern global access method and strict wave access.
#pragma version=7.0
/////////////////////////////////////////////
// Spectromicsroscopy beamline Elettra
// Alexei Barinov
// Marek Kopciuszynski marek.kopciuszynski@elettra.eu
/////////////////////////////////////////////

// SPEM init file
// Load the procedures and setup the environment
#include "SPEM_Interface" // Procedures responsible for user interface
#include "SPEM_LoadData"  // Loading the data
#include "SPEM_PlotData"  // Prepare plots
#include "SPEM_ProcessData" // Process the data
#include "SPEM_Kspace"		// Transform maps to k space
#include  <All IP Procedures> // Link to standard image processing tool (Igor)
///////////////////////////////
// String constants
///////////////////////////////
//Script data folder
Strconstant ksSPMDataF="root:SPMData:"
Strconstant ksDataFolder="root:"
//////////////////////////////////////
//Wave with attributes for loaded data
StrConstant ksTextAttributes="root:SPMData:txtAttributes"
//Titles and utility names of the HDF-SPM attributes
Strconstant ksAtts="Dataset name;Sample ID;Acquisition Type;Name 0;Units 0;Offset 0;Delta 0;Name 1;Units 1;Offset 1;Delta 1;Name 2;Units 2;Offset 2;Delta 2;X;Y;Z;R;aT;aP;Det;MCP;Lens Mode;PE;Dwell;Scans;Initial channel;Final channel;Initial channel Y;Final channel Y;T;DT;P;Ring Energy;Gap;Photon Energy;Ring current, start;Ring Current, end;Date time start;Date Time End;Calibration;defocus;Name3;Units3;Offset3;Del3;"
//Every attribute name expressed as a constant
StrConstant ksSampleID="Sample ID", ksAcqType="Acquisition Type",ksDatasetID="Dataset name",ksDetector="Det",ksClbW="Calibration"
StrConstant ksStarted="Date time start",ksFinished="Date Time End"
StrConstant ksT="T",ksDT="DT",ksX="X",ksY="Y",ksZ="Z",ksR="R",ksP="P",ksd="defocus"
StrConstant ksName0="Name 0",ksUnits0="Units 0",ksOffset0="Offset 0",ksDelta0="Delta 0",ksName1="Name 1",ksUnits1="Units 1",ksOffset1="Offset 1",ksDelta1="Delta 1",ksName2="Name 2",ksUnits2="Units 2",ksOffset2="Offset 2",ksDelta2="Delta 2",ksName3="Name3",ksUnits3="Units3",ksOffset3="Offset3",ksDelta3="Del3"
StrConstant ksaT="aT",ksaP="aP"
StrConstant ksInitialCh="Initial channel", ksFinalCh="Final channel", ksScans="Scans",ksKEmin="Initial channel Y",ksKEmax="Final channel Y",ksPE="PE",ksDwell="Dwell",kcMCP="MCP",ksLensMode="Lens Mode"
StrConstant ksRingE="Ring Energy",ksRingCur0="Ring current, start",ksRingCur1="Ring Current, end",ksGap="Gap",ksPhoton="Photon Energy"
//Acquisition types
StrConstant ksImage="Image",ksPolar="Polar Scan",ksFocus="Focus Scan",ksSpectrum="Spectrum",ksEnergy="Energy Map",ksClb="Calibration", ksImage4d="4D", ksAngularScan = "Angular Scan"


/////////////////////////////////////
// Set the default color map for ARPES data
StrConstant ksColorMap = "viridis"


// Setup the environment
proc SPEM_Init()
	SetDataFolder ksDataFolder

	NewDataFolder /O/S root:SPMData

	variable/G AreaSize = 1
	variable/G SpectraIndex = 0
	variable/G ActiveClb = -1

	// Default values for Fermi Energy and P/T values for normal emission
	variable/g gef74 = 70.8
	variable/g gef27 = 23.7
	variable/g gT0 = 0
	variable/g gP90 = 88
	variable/g gPgon = 0

	SetDataFolder ksDataFolder
	// Load color map from matplotlib tables
	LoadWave/O/Q/P = IgorUserFiles ":User Procedures:ElettraSPEM:viridis.ibw"

	print " "
	print "==========================================="
	print " Elettra SPEM routines loaded successfully "
	print "==========================================="
	print " "
	print "The latest version and instructions could be found here:"
	print "https://github.com/mkopciuszynski/ElettraSPEM"
	print " "
end
