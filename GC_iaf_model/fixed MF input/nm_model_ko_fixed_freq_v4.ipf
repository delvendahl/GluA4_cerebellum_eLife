#pragma rtGlobals=3		// Use modern global access method and strict wave access.

// Simulation of GC with fixed MF input
//
// Definition of model parameters:
//
// BEGIN definitions
static constant TIMESTEP = 0.05	// time step (ms)
static constant STIMSTART = 50		// start of stimulation (ms)
static constant STIMNO = 10			// number of stimuli
static constant DURATION = 200		// simulation time (ms)
static constant BINOMIAL_N = 10	// N for binomial ms)
static constant SIMULATIONS = 20	// number of simulations
static constant DIAMETER = 11.7	// model diameter (µm), gives Cm=4.3005pF
static constant GTONIC = 0.0		// tonic GABA conductance (nS)
static constant GLEAK = 0.333		// leak conductance (nS)
static constant EGABA = -65 		// GABA reversal potential (mV)
static constant ELEAK = -100		// leak reversal potential (mV)
// END definitions


function init()
	NMTabAdd( "Pulse", "" )
	NMTabRemove( "Fit" )
	NMTabAdd( "Model", "" )
	NMTabRemove( "Event" )
	setrandomseed 0.4
end

Function run()
	Init()
	
	variable numV = 4
	Make/O/N=(numV) voltages = -100 + p*10
	variable numFreq = 3
	variable increment = 100
	Make/O/N=(numFreq) MF_Freq = increment + p * increment
	
	Make/O/D/N=(numFreq,numV) GC_Freq, GC_Freq_SD, Fano, VRE, Delay, Delay_SD

	variable i,j
	for(i = 0; i < numV; i += 1)
		Printf "Model set at %+.0f mV\r", voltages[i]
		Set_model(voltages[i])
		
		for(j = 0; j < numFreq; j += 1)
			Printf "Running model with: %.1f Hz MF stimulation\r", MF_Freq[j]
	
			wave results = run_model(MF_Freq[j], voltages[i])
			GC_Freq[j][i] = results[0]
			GC_Freq_SD[j][i] = results[1]
			Fano[j][i] = results[2]
			Delay[j][i] = results[3]
			Delay_SD[j][i] = results[4]
	
			string PrefStr
			sprintf PrefStr "V%03.0f_F%04.0f_", abs(voltages[i]),MF_Freq[j]
			VRE[j][i] = GetVanRossumError(PrefStr)
			printf "\r\r"
		endfor	
	endfor	
	KillWaves/Z results	
	edit MF_Freq, GC_Freq, GC_Freq_SD, Delay, Delay_SD, Fano, VRE
	
	//clean up
	string wavestokill = wavelist("Sim_*", ";", "")
	wavestokill += wavelist("Model_*", ";", "")
	wavestokill += wavelist("gGABA*", ";", "")
	for(i = 0; i < itemsinlist(wavestokill); i += 1)
		wave w = $stringfromlist(i, wavestokill)
		Killwaves/Z w
	endfor
	
	//plot results
	Plot2d_Wave(GC_Freq, MF_Freq, "Firing (Hz)")
	Plot2d_Wave(Delay, MF_Freq, "Delay (ms)")

	
end

Static function Plot2d_Wave(ywave,xwave,ylabel)
	wave ywave,xwave
	string ylabel
	
	if(dimsize(ywave,1) == 1)
		Display ywave[][0] vs xwave
	else
		Display ywave[][1] vs xwave
		AppendToGraph ywave[][2] vs xwave
		AppendToGraph ywave[][3] vs xwave
	endif
	ModifyGraph mode=4,marker=19,nticks=2
	ModifyGraph width=70.8661,height=70.8661,gfSize=7,expand=2
	SetAxis/A/N=1/E=1 left
	Label bottom "MF stimulation (Hz)"
	Label left ylabel
	string traces = TraceNameList("",";",1)
	ModifyGraph rgb($stringfromlist(0,traces))=(48059,48059,48059)
	ModifyGraph rgb($stringfromlist(1,traces))=(30583,30583,30583)	
	ModifyGraph rgb($stringfromlist(2,traces))=(0,0,0)

End

Static Function/WAVE run_model(frequency, voltage)
	variable frequency,voltage
	
	NMSet( tab="Pulse" )
	NM_make_Conduct_Waves_FixedFreq(frequency)

	NMSet( tab="Model" )
	NMModelRun()
	
	NMSet( tab="Spike" )
	Variable startTime = STIMSTART
	Variable endTime = STIMSTART + STIMNO*(1/frequency)*1000 + 10	// added 10 ms
	NMSpikeSet( xbgn=startTime )
	NMSpikeSet( xend=endTime )
	NMSpikeRasterComputeAll( chanSelectList="A", waveSelectList="All", displayMode=1, delay=0 )
	NMSpikeRate( folder="Spike_Sim_Vmem_All_A", xRasterList="SP_RX_SimVmemAll_A0", xbgn=startTime, xend=endTime )
	
	Wavestats/Q root:nmFolder0:Spike_Sim_Vmem_All_A:Rate_SP_RX_SimVmemAll_A0
	variable rate = V_avg * 1000
	variable rate_sd = V_sdev * 1000
	
	Wave SpikeDelays = CalcSpikeDelay()
	Wavestats/Q SpikeDelays
	variable delay = V_avg
	variable delay_sd = V_sdev
	Duplicate/O SpikeDelays, $("SpikeDelays_" + num2str(abs(voltage)) + "mV_" +  num2str(frequency) + "Hz")
	
	Make/O/D/FREE/N=5 wOut = {rate, rate_sd, (rate_sd^2) / rate, delay, delay_sd}
	
	Printf "Simulation run with MF input at %.0f Hz\r",frequency
	Printf "GC frequency: %.5f ± %.5f Hz\r", rate, rate_sd
	Printf "First AP delay: %.5f ± %.5f ms\r", delay, delay_sd
	Printf "Fano factor: %.5f Hz\r",(rate_sd^2)/rate 
	
	NMSet( tab="Main" )
	
	string PrefStr
	sprintf PrefStr "V%03.0f_F%04.0f_", abs(voltage),frequency
	NMMainDuplicate( newPrefix=PrefStr )
	sort_into_folder(PrefStr)
	
	return wOut
End

Static Function Set_model(vmem)
	variable vmem
	
	variable surface = DIAMETER^2*pi
	variable gTonicDensity = GTONIC / surface
	variable gLeakDensity = GLEAK / surface
	variable iClamp = GTONIC*(vmem-EGABA) + GLEAK*(vmem-ELEAK)
	
	NMSet( tab="Model" )
	NMModelSet( modelSelect="IAF_AdEx" )
	NMModelVarSet( "SimulationTime", DURATION )
	NMModelVarSet( "TimeStep", TIMESTEP )
	NMModelVarSet( "Temperature", 25 )
	NMModelVarSet( "V0", vmem )
	NMModelVarSet( "Diameter", DIAMETER )
	NMModelVarSet( "gTonicGABADensity", gTonicDensity )
	NMModelStrSet( "gGABA_WavePrefix", "gGABA__" )
	NMModelVarSet( "eGABA", EGABA )
	NMModelStrSet( "gNMDA_BlockFxnList", "GC_Schwartz2012" )
	NMModelVarSet( "gLeakDensity", gLeakDensity )
	NMModelVarSet( "eLeak", ELEAK )
	NMModelVarSet( "AP_Threshold", -45 )
	NMModelVarSet( "AP_Peak", 32 )
	NMModelVarSet( "AP_Reset", -75 )
	NMModelVarSet( "AP_Refrac", 1.5 )
	NMModelVarSet( "W_Tau", 75 )
	NMModelVarSet( "W_A", 0 )
	NMModelVarSet( "W_B", 0.35 )
	
	NMModelVarSet( "iClampAmp", iClamp )
	NMModelVarSet( "iClampOnset", 0 )
	NMModelVarSet( "iClampDuration", DURATION )
	
	NMModelVarSet( "LastSimulation", SIMULATIONS-1 )
	
	NMModelTabUpdate()
End
	

Static Function NM_make_Conduct_Waves_FixedFreq(Freq) // code for generating MF-GC conductance waves
	Variable Freq // Hz, freq of MF input
	
	Freq/= 1000 //kHz
	Variable interval = 1/Freq	// ISI in ms
		
	Variable SavePlasticityWavesAMPA = 0 // ( 0 ) no ( 1 ) yes
	Variable SavePlasticityWavesNMDA = 0
	
	Variable tbgn = STIMSTART
	Variable tend = STIMSTART + interval*(STIMNO-1)
	
	String gAMPA_Prefix = "gAMPA_"
	String gGABA_Prefix = "gGABA_"
	String gNMDA_Prefix = "gNMDA_"
	
	String directList1 = KO_NMPulseExp_GC_AMPAdirect( select = 0 ) + KO_NMPulseTrainRP_GC_AMPAdirect()
	String directList2 = KO_NMPulseExp_GC_AMPAdirect( select = 1 ) + KO_NMPulseTrainRP_GC_AMPAdirect()
	String spillList1 = KO_NMPulseExp_GC_AMPAspill( select = 0 ) + KO_NMPulseTrainRP_GC_AMPAspill()
	String spillList2 = KO_NMPulseExp_GC_AMPAspill( select = 1 ) + KO_NMPulseTrainRP_GC_AMPAspill()
	String spillList3 = KO_NMPulseExp_GC_AMPAspill( select = 2 ) + KO_NMPulseTrainRP_GC_AMPAspill()
	String nmdaList1 = KO_NMPulseExp_GC_NMDA( select = 0 ) + KO_NMPulseTrainRP_GC_NMDA()
	String nmdaList2 = KO_NMPulseExp_GC_NMDA( select = 1 ) + KO_NMPulseTrainRP_GC_NMDA()
		
	NMPulseSet( wavePrefix=gAMPA_Prefix, update=1 )
	
	NMPulseSet( xunits="msec" )
	NMPulseSet( yunits="nS" )
	NMPulseSet( numWaves=SIMULATIONS )
	NMPulseSet( dx=TIMESTEP )
	NMPulseSet( waveLength=DURATION )
	NMPulseConfigRemove( all=1 )
	
	NMPulseSet(  wavePrefix=gNMDA_Prefix, update=1 ) // config # 2
	NMPulseSet(  wavePrefix=gNMDA_Prefix, update=1 )
	NMPulseSet( xunits="msec" )
	NMPulseSet( yunits="nS" )
	NMPulseSet( numWaves=SIMULATIONS )
	NMPulseSet( dx=TIMESTEP )
	NMPulseSet( waveLength=DURATION )
	NMPulseConfigRemove( all=1 )
	
	// train parameters
	String paramList = "train=fixed" + ";tbgn=" + num2str( tbgn ) + ";tend=" + num2str( tend ) + ";interval=" + num2str( interval )
	
	// gAMPA (direct + spillover)	
	NMPulseSet( wavePrefix=gAMPA_Prefix, update=1 )
	NMPulseSet( wavePrefix=gAMPA_Prefix, update=1 )
	
	NMPulseConfigAdd( "wave=all" + ";" + paramList + ";" + directList1 + ";binomialN=" + num2str(BINOMIAL_N) + ";" )
	NMPulseConfigAdd( "wave=all" + ";" + paramList + ";" + directList2 + ";binomialN=" + num2str(BINOMIAL_N) + ";" )
	NMPulseConfigAdd( "wave=all" + ";" + paramList + ";" + spillList1 + ";binomialN=" + num2str(BINOMIAL_N) + ";" )
	NMPulseConfigAdd( "wave=all" + ";" + paramList + ";" + spillList2 + ";binomialN=" + num2str(BINOMIAL_N) + ";" )
	NMPulseConfigAdd( "wave=all" + ";" + paramList + ";" + spillList3 + ";binomialN=" + num2str(BINOMIAL_N) + ";" )
	
	// gNMDA
	NMPulseSet(  wavePrefix=gNMDA_Prefix )
	NMPulseSet(  wavePrefix=gNMDA_Prefix )
		
	NMPulseConfigAdd( "wave=all" + ";" + paramList + ";" + nmdaList1 + ";binomialN=" + num2str(BINOMIAL_N) + ";" )
	NMPulseConfigAdd( "wave=all" + ";" + paramList + ";" + nmdaList2 + ";binomialN=" + num2str(BINOMIAL_N) + ";" )
	
	NMPulseSet( wavePrefix=gAMPA_Prefix, update=1 )
	NMPulseSet( wavePrefix=gAMPA_Prefix, update=1 )
	NMConfigVarSet( "Pulse", "SavePlasticityWaves", SavePlasticityWavesAMPA )
	NMPulseExecute()

	NMPulseSet( wavePrefix=gNMDA_Prefix, update=1 )
	NMPulseSet( wavePrefix=gNMDA_Prefix, update=1 )
	NMConfigVarSet( "Pulse", "SavePlasticityWaves", SavePlasticityWavesNMDA )
	NMPulseExecute()

	// gGABA (not used)
	Variable i
	for ( i = 0 ; i < SIMULATIONS; i += 1 )
		Variable npnts = numpnts( $gAMPA_Prefix + num2istr( i ) )
		String wName = gGABA_Prefix + num2istr( i )
	
		Make/O/N= ( npnts) $wName = GTONIC
		SetScale/P x 0, TIMESTEP, "", $wName
	endfor
End


// *************************************
// STATIC functions generating parameter lists for AMPA & NMDA

STATIC Function /S WT_NMPulseExp_GC_AMPAdirect( [ select, p ] ) // Billings 2014
	Variable select // 0 or 1
	STRUCT NMPulseExp &p
	STRUCT NMPulseExp pp
	
	Variable amp = NaN
	
	if ( !ParamIsDefault( p ) )
		pp = p
	endif
	
	pp.amp1 = -1	
	pp.tau1 = 0.3274 // ms
	pp.amp2 = 1
	
	switch( select )
		case 0:
			amp = 1.4/BINOMIAL_N // nS
			pp.tau2 =  0.5 // ms
			break
		case 1:
			amp = 0.55/BINOMIAL_N // nS
			pp.tau2 = 4 // ms
			break
	endswitch
	
	pp.amp3 = 0
	pp.amp4 = 0
	
	return "pulse=exp;amp=" + num2str( amp ) + ";" + NMPulseExpParamList( pp )
	
End // WT_NMPulseExp_GC_AMPAdirect

STATIC Function /S WT_NMPulseTrainRP_GC_AMPAdirect( [ p ] )
	STRUCT NMPulseTrainRP &p
	STRUCT NMPulseTrainRP pp
	
	if ( !ParamIsDefault( p ) )
		pp = p
	endif
	
	// Billings 2014
	pp.Rinf = 1
	pp.Rmin = 0
	pp.tauR = 85
	pp.Pinf = 0.45
	pp.Pmax = inf
	pp.tauP = -1 // this turns P off
	pp.Pscale = 0
	
	return NMPulseTrainRPparamList( pp )

End // WT_NMPulseTrainRP_GC_AMPAdirect

STATIC Function /S WT_NMPulseExp_GC_AMPAspill( [ select, p ] ) // Billings 2014
	Variable select // 0 or 1 or 2
	STRUCT NMPulseExp &p
	STRUCT NMPulseExp pp
	
	Variable amp = NaN
	
	if ( !ParamIsDefault( p ) )
		pp = p
	endif
	
	pp.amp1 = -1
	pp.tau1 = 0.9 // ms
	pp.amp2 = 1
	
	switch( select )
		case 0:
			amp = 0.1/BINOMIAL_N // nS
			pp.tau2 = 0.8 // ms
			break		
		case 1:
			amp = 0.25/BINOMIAL_N // nS
			pp.tau2 = 4 // ms
			break	
		case 2:
			amp = 0.15/BINOMIAL_N // nS
			pp.tau2 = 30 // ms
			break
		default:
			return ""			
	endswitch
	
	return "pulse=exp;amp=" + num2str( amp ) + ";" + NMPulseExpParamList( pp )
	
End // WT_NMPulseExp_GC_AMPAspill

 STATIC Function /S WT_NMPulseTrainRP_GC_AMPAspill( [ p ] )
	STRUCT NMPulseTrainRP &p
	STRUCT NMPulseTrainRP pp
	
	if ( !ParamIsDefault( p ) )
		pp = p
	endif
	
	// Billings 2014
	pp.Rinf = 1
	pp.Rmin = 0
	pp.tauR = 25
	pp.Pinf = 0.45
	pp.Pmax = inf
	pp.tauP = -1 // this turns P off
	pp.Pscale = 0
	
	return NMPulseTrainRPparamList( pp )

End // WT_NMPulseTrainRP_GC_AMPAspill

STATIC Function /S WT_NMPulseExp_GC_NMDA( [ select, p ] ) // Billings 2014
	Variable select // 0 or 1
	STRUCT NMPulseExp &p
	STRUCT NMPulseExp pp
	
	Variable amp = NaN
	
	if ( !ParamIsDefault( p ) )
		pp = p
	endif
	
	pp.amp1 = -1
	pp.tau1 = 2.5 // ms
	pp.amp2 = 1
	
	switch( select )
		case 0:
			amp = 0.3/BINOMIAL_N // nS
			pp.tau2 = 30 // ms
			break
		case 1:
			amp = 0.14/BINOMIAL_N // nS
			pp.tau2 = 70 // ms
			break
		default:
			return ""
	endswitch
	
	return "pulse=exp;amp=" + num2str( amp ) + ";" + NMPulseExpParamList( pp )
	
End // WT_NMPulseExp_GC_NMDA

STATIC Function /S WT_NMPulseTrainRP_GC_NMDA( [ p ] )
	STRUCT NMPulseTrainRP &p
	STRUCT NMPulseTrainRP pp
	
	if ( !ParamIsDefault( p ) )
		pp = p
	endif
	
	// Billings 2014
	pp.Rinf = 1
	pp.Rmin = 0
	pp.tauR = 50
	pp.Pinf = 0.45
	pp.Pmax = inf
	pp.tauP = 15
	pp.Pscale = 0.2
	
	return NMPulseTrainRPparamList( pp )

End // WT_NMPulseTrainRP_GC_NMDA

///////////////////////

STATIC Function /S KO_NMPulseExp_GC_AMPAdirect( [ select, p ] ) // Billings 2014
	Variable select // 0 or 1
	STRUCT NMPulseExp &p
	STRUCT NMPulseExp pp

	Variable amp = NaN

	if ( !ParamIsDefault( p ) )
		pp = p
	endif

	pp.amp1 = -1
	pp.tau1 = 0.45 // ms
	pp.amp2 = 1
	pp.amp3 = 0
	pp.amp4 = 0

	switch( select )
		case 0:
			amp = 0.196/BINOMIAL_N // nS
			pp.tau2 =  0.9 // ms
			break
		case 1:
			amp = 0.077/BINOMIAL_N // nS
			pp.tau2 = 7.2 // ms
			break
	endswitch

	return "pulse=exp;amp=" + num2str( amp ) + ";" + NMPulseExpParamList( pp )
End // KO_NMPulseExp_GC_AMPAdirect

STATIC Function /S KO_NMPulseTrainRP_GC_AMPAdirect( [ p ] )
	STRUCT NMPulseTrainRP &p
	STRUCT NMPulseTrainRP pp

	if ( !ParamIsDefault( p ) )
		pp = p
	endif

	pp.Rinf = 1
	pp.Rmin = 0
	pp.tauR = 57
	pp.Pinf = 0.4
	pp.Pmax = inf
	pp.tauP = -1 // this turns P off
	pp.Pscale = 0

	return NMPulseTrainRPparamList( pp )
End // KO_NMPulseTrainRP_GC_AMPAdirect

STATIC Function /S KO_NMPulseExp_GC_AMPAspill( [ select, p ] ) // Billings 2014
	Variable select // 0 or 1 or 2
	STRUCT NMPulseExp &p
	STRUCT NMPulseExp pp

	Variable amp = NaN

	if ( !ParamIsDefault( p ) )
		pp = p
	endif

	pp.amp1 = -1
	pp.tau1 = 1.1 // ms
	pp.amp2 = 1

	switch( select )
		case 0:
			amp = 0.014/BINOMIAL_N // nS
			pp.tau2 = 1.44 // ms
			break
		case 1:
			amp = 0.035/BINOMIAL_N // nS
			pp.tau2 = 7.2 // ms
			break
		case 2:
			amp = 0.021/BINOMIAL_N // nS
			pp.tau2 = 54 // ms
			break
		default:
			return ""
	endswitch

	return "pulse=exp;amp=" + num2str( amp ) + ";" + NMPulseExpParamList( pp )
End // KO_NMPulseExp_GC_AMPAspill

STATIC Function /S KO_NMPulseTrainRP_GC_AMPAspill( [ p ] )
	STRUCT NMPulseTrainRP &p
	STRUCT NMPulseTrainRP pp

	if ( !ParamIsDefault( p ) )
		pp = p
	endif

	pp.Rinf = 1
	pp.Rmin = 0
	pp.tauR = 17
	pp.Pinf = 0.4
	pp.Pmax = inf
	pp.tauP = -1 // this turns P off
	pp.Pscale = 0

	return NMPulseTrainRPparamList( pp )
End // KO_NMPulseTrainRP_GC_AMPAspill

STATIC Function /S KO_NMPulseExp_GC_NMDA( [ select, p ] ) // Billings 2014
	Variable select // 0 or 1
	STRUCT NMPulseExp &p
	STRUCT NMPulseExp pp

	Variable amp = NaN

	if ( !ParamIsDefault( p ) )
		pp = p
	endif

	pp.amp1 = -1
	pp.tau1 = 2.5 // ms
	pp.amp2 = 1

	switch( select )
		case 0:
			amp = 0.45/BINOMIAL_N // nS
			pp.tau2 = 30 // ms
			break
		case 1:
			amp = 0.21/BINOMIAL_N // nS
			pp.tau2 = 70 // ms
			break
		default:
			return ""
	endswitch

	return "pulse=exp;amp=" + num2str( amp ) + ";" + NMPulseExpParamList( pp )
End // KO_NMPulseExp_GC_NMDA

STATIC Function /S KO_NMPulseTrainRP_GC_NMDA( [ p ] )
	STRUCT NMPulseTrainRP &p
	STRUCT NMPulseTrainRP pp

	if ( !ParamIsDefault( p ) )
		pp = p
	endif

	pp.Rinf = 1
	pp.Rmin = 0
	pp.tauR = 50
	pp.Pinf = 0.45
	pp.Pmax = inf
	pp.tauP = 15
	pp.Pscale = 0.2

	return NMPulseTrainRPparamList( pp )
End // KO_NMPulseTrainRP_GC_NMDA

//////////////////////////////////

Static function GetVanRossumError(folderStr)
	string folderStr

	DFREF savedfr = GetDataFolderDFR()
	SetDataFolder root:nmFolder0:$(folderstr):

	string rWaves = WaveList(folderStr + "Sim*", ";", "")
	variable numWaves = itemsinlist(rWaves)

	Make/D/O/N=(numWaves) $(folderStr + "van_rossum") /WAVE=vre
	variable dummy,freq
	sscanf folderStr, "V%d_F%d", dummy,freq

	variable i
	for(i = 0; i < numWaves; i += 1)
		wave w = $stringfromlist(i, rWaves)
		Make/O/D/N=0 $("vRE_GC"+num2istr(i)) /WAVE=gc
		Findlevels/Q/D=gc/M=0.5/R=(STIMSTART,STIMSTART + STIMNO*(1/freq*1000)) w, -10

		Make/O/D/N=(STIMNO) $(folderStr + "train") /WAVE=mf
		mf = STIMSTART + p*(1/freq*1000)
		
		wave target = wave_from_spiketimes(mf)
		wave result = wave_from_spiketimes(gc)

		vre[i] = van_rossum(target, result, 10)
	endfor

	string kList = wavelist("vRE_*", ";", "")
	for(i = 0; i < itemsinlist(kList); i += 1)
		KillWaves/Z $stringfromlist(i, kList)
	endfor

	Wavestats/Q vre
	Printf "van Rossum error: %.5f ± %.5f\r", v_avg, v_sdev

	SetDataFolder savedfr
	
	return v_avg
end

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

static function/WAVE wave_from_spiketimes(inwave)
	wave inwave

	Make/O/N=(DURATION * 200) $(nameofwave(inwave) + "_spikes") / wave = out
	setscale/P x 0, (1 / 200), "ms", out

	variable i
	for(i = 0; i < numpnts(inwave); i += 1)
		variable pnt = x2pnt(out, inwave[i])
		out[pnt] = 1
	endfor

	return out
end

Static Function sort_into_folder(PrefStr)
	string PrefStr

	DFREF saveDFR = GetDataFolderDFR()
	
	NewDataFolder/O/S $PrefStr
	DFREF newDFR = GetDataFolderDFR()
	SetDataFolder saveDFR
	
	string wList = WaveList("gAMPA*", ";", "")
	wList += WaveList("gNMDA*", ";", "")
	wList += WaveList("SpikeDelays_*", ";", "")
	wList += WaveList(PrefStr+"*", ";", "")
	variable i
	for(i = 0; i < itemsinlist(wList); i += 1)
		wave w = $stringfromlist(i, wList)
		Movewave w, newDFR
	endfor
End

Function/WAVE CalcSpikeDelay()
	wave xdat = $("root:nmFolder0:Spike_Sim_Vmem_All_A:SP_RX_SimVmemAll_A0")
	wave ydat = $("root:nmFolder0:Spike_Sim_Vmem_All_A:SP_RY_SimVmemAll_A0")

	variable numcells = wavemax(ydat)+1
	Make/D/FREE/N=(numcells) SpikeDelays
	
	variable i, j = 0
	variable cellno = 0
	for(i=0; i < numpnts(ydat); i += 1)	
		if(ydat[i] == cellno)
			SpikeDelays[j] = xdat[i] - STIMSTART
			j += 1
			cellno += 1
		endif
	endfor
	
	return SpikeDelays
End