#pragma rtGlobals=3

// Load SPEM procedures 
Menu "Data"
	"-"
	"Load SPEM", Execute/P "INSERTINCLUDE \"SPEM_Init\""; Execute/P "COMPILEPROCEDURES "; Execute/P "SPEM_Init()"
End



