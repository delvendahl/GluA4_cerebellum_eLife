#pragma rtGlobals=3		// Use modern global access method and strict wave access.

// calculates the GC threshold (=number of MFs required to drive spiking)
// presynaptic spikes are integrated over 25 ms before a given GC spike (variable t_end)

function get_threshold()

	Wave Runs
	variable numFreq = numpnts(Runs)

	DFREF savedfr = GetDataFolderDFR()
	string wList = "", mfList = ""
	
	variable t_start = 0
	variable t_end = 25
	
	Make/O/D/N=(50,81) active_mfs
	Make/O/D/N=(50,81) gc_spikes
	
	variable i
	
	for(i = 0; i < numFreq; i += 1)
		string folderStr
		sprintf folderStr "P%03.0f", Runs[i]-1
		
		SetDataFolder root:nmFolder0:$(folderstr):
		string rWaves = WaveList("Sim_*", ";", "")
		variable numWaves = itemsinlist(rWaves)
		
		variable j
		for(j = 0; j < numWaves; j += 1)
			wave w = $stringfromlist(j, rWaves)
			Findlevels/Q/M=0.5/R=(t_start,t_end) w, -10
			gc_spikes[i][j] = V_LevelsFound
			
			wave conn =  root:nmFolder0:GC_to_MF
			variable MFs = 0
			variable k
			for(k = 0; k<4; k+=1)
				wave mf_wave = $("train" + num2str(conn[j][k]) + "_" + num2str(i))		
				Extract/FREE mf_wave, temp, (mf_wave > t_start && mf_wave< t_end)
				MFs += (dimsize(temp,0)) ? 1 : 0
			endfor
			active_mfs[i][j] = MFs

		endfor 
		
		SetDataFolder savedfr
	endfor

	redimension /N=(dimsize(active_mfs,0)*dimsize(active_mfs,1)) active_mfs
	redimension /N=(dimsize(gc_spikes,0)*dimsize(gc_spikes,1)) gc_spikes
	
	extract/O active_mfs, thresholds, (gc_spikes >= 1)
	
	Make/N=6/O thresholds_Hist;DelayUpdate
	Histogram/P/C/B={0,1,6} thresholds,thresholds_Hist
	
	Display thresholds_Hist
	ModifyGraph mode=5,offset={-0.5,0},hbFill=4
	ModifyGraph offset(thresholds_Hist)={-1,0}
	
	Make/N=6/O thresholds_cumHist;DelayUpdate
	Histogram/CUM/P/B={0,1,6} thresholds,thresholds_cumHist;DelayUpdate
	AppendToGraph /R thresholds_cumHist
	ModifyGraph offset(thresholds_cumHist)={-0.5,0}
	
end