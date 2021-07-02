#pragma rtGlobals=3		// Use modern global access method and strict wave access.

// functions to export spike time data to .txt
// data are saved to folder where .pxp file is located 
// run "export_all()" to iterate over folders

// export spike times 
function ExportSpikeTimes(folderStr)
	string folderStr

	string prefix = "WT_"
	DFREF savedfr = GetDataFolderDFR()

	SetDataFolder root:nmFolder0:$(folderstr):

	string rWaves = WaveList(folderStr + "Sim*", ";", "")
	variable numWaves = itemsinlist(rWaves)

	PathInfo home
	NewPath/O/Q dest, S_path
	string filenamestring = prefix + folderstr + ".txt"
	
	variable i
	for(i = 0; i < numWaves; i += 1)
		wave w = $stringfromlist(i, rWaves)
		Make/O/D/N=0 $("temp_GC"+num2istr(i)) /WAVE=gc
		Findlevels/Q/D=gc/M=0.5 w, -10

		string tList = WaveList("train" + num2istr(i) + "_*", ";", "")
		concatenate/NP/O tList, $("temp_MF"+num2istr(i)) /WAVE=mf
		sort mf, mf
		
		makeUnique(mf)
		
		redimension/N=(1,numpnts(mf)) mf
		redimension/N=(1,numpnts(gc)) gc
		
		save/O/A=2/J/M="\r\n"/P=dest mf as filenamestring
		save/O/A=2/J/M="\r\n"/P=dest gc as filenamestring
		
	endfor 

	string kList = wavelist("temp_*", ";", "")
	for(i = 0; i < itemsinlist(kList); i += 1)
		KillWaves/Z $stringfromlist(i, kList)
	endfor

	SetDataFolder savedfr
end

// wrapper function to export data for all frequencies
function export_all()

	wave MF_Freq
	variable numFreq = numpnts(MF_Freq)
	variable i
	
	for(i = 0; i < numFreq;  i += 1)
		string folderStr
		sprintf folderStr "P%04.0f_", MF_Freq[i]
		
		ExportSpikeTimes(folderStr)
	endfor
end

// remove duplicate spike times
static function makeUnique(inWave)
	wave inWave
	
	variable num = numpnts(inWave)
	variable i
	for(i = 0; i < (num - 1); i += 1)
		inWave[i] = (inWave[i] == inWave[i+1]) ? NaN : inWave[i]
	endfor
	Wavetransform ZapNaNs inWave
end	