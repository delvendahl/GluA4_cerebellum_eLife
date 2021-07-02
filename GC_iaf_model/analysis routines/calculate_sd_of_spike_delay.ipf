#pragma rtGlobals=3		// Use modern global access method and strict wave access.

// calculates standard deviation of first AP delay
// loops through folders created by the simulation 
// naming convention: "Pxxxx_" with xxxx being stimulation frequency

function calculate_delay_stdev()

	wave Delay_SD
	dfref dfSav = GetDataFolderDFR()
	variable i
	for(i = 0; i < numpnts( Delay_SD ); i += 1)
		string PrefStr
		sprintf PrefStr "P%04.0f_", 20+i*20
		setdatafolder $PrefStr
		wave w = $("SpikeDelays_"+ num2str(20+i*20) + "Hz")
		wavestats/Q w
		Delay_SD[i] = v_sdev
		SetDataFolder dfSav
	endfor

end
