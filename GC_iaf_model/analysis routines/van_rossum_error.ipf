#pragma rtGlobals=3		// Use modern global access method and strict wave access.

// functions to calculate the van Rossum error between presynaptic and postsynaptic spike trains
// uses as exponential kernel (default: 10ms)
// result is printed to history and copied to clipboard
// run as "GetVanRossumError(folderStr)"

// main function to calculate van rossum distance
// params: folderStr (string), name of data folder to operate on, must contain simulation data
function GetVanRossumError(folderStr)
	string folderStr

	variable kernel_tau = 10

	DFREF savedfr = GetDataFolderDFR()

	SetDataFolder root:nmFolder0:$(folderstr):

	string rWaves = WaveList(folderStr + "Sim*", ";", "")
	variable numWaves = itemsinlist(rWaves)

	Make/D/O/N=(numWaves) $(folderStr + "van_rossum") /WAVE=vre

	variable i
	for(i = 0; i < numWaves; i += 1)
		wave w = $stringfromlist(i, rWaves)
		Make/O/D/N=0 $("vRE_GC"+num2istr(i)) /WAVE=gc
		Findlevels/Q/D=gc/M=0.5/R=(200,1200) w, -10

		string tList = WaveList("train" + num2istr(i) + "_*", ";", "")
		concatenate/NP/O tList, $("vRE_MF"+num2istr(i)) /WAVE=mf
		sort mf, mf

		wave target = wave_from_spiketimes(mf)
		wave result = wave_from_spiketimes(gc)

		vre[i] = van_rossum(target, result, kernel_tau)
	endfor

	string kList = wavelist("vRE_*", ";", "")
	for(i = 0; i < itemsinlist(kList); i += 1)
		KillWaves/Z $stringfromlist(i, kList)
	endfor

	Wavestats/Q vre

	Print folderStr
	Printf "van Rossum error: %.5f � %.5f\r\r", v_avg, v_sdev

	string clip
	sprintf clip "%.8f\t%.8f", v_avg, v_sdev
	putscraptext clip

	SetDataFolder savedfr
end


// create a continuous wave given spiek time data
// params: inwave, wave containing spike timings
// returns continuous data wave, length 1500 ms

static function/WAVE wave_from_spiketimes(inwave)
	wave inwave

	variable length = 1500	// length in ms

	Make/O/N=(length * 200) $(nameofwave(inwave) + "_spikes") / wave = out
	setscale/P x 0, (1 / 200), "ms", out

	variable i
	for(i = 0; i < numpnts(inwave); i += 1)
		variable pnt = x2pnt(out, inwave[i])
		out[pnt] = 1
	endfor

	return out
end

// calculate the van rossum distance between two waves
// params: wave1, wave2 are input waves; kernelTau is the exponential kernel to use
// returns van rossum distance (scalar)
static function van_rossum(wave1, wave2, kernelTau)
	wave wave1,wave2	// scaled in ms
	variable kernelTau	// ms

	wave a = convolve_wave(wave1, kernelTau)
	wave b = convolve_wave(wave2, kernelTau)

	duplicate/O/FREE b difference
	difference -= a

	Duplicate/O/FREE difference temp
	temp = difference^2

	variable vanrossum = (1/kernelTau) * area(temp)

	return  sqrt(vanrossum)
end

// convolve a wave with a given exponential kernel
// params: inwave is the input wave; kernelTau is the exponential kernel to use
// returns convolved wave
static function/WAVE convolve_wave(inwave, kernelTau)
	wave inwave
	variable kernelTau

	Make/O/FREE/N=(numpnts(inwave)/10) kernel
	setscale/P x 0, (1 / 200), "ms", kernel
	kernel = exp(- x / kernelTau)

	duplicate/O inwave $(nameofwave(inwave) + "_conv") /WAVE=w
	convolve kernel, w

	deletepoints numpnts(inwave), inf, w

	return w
end
