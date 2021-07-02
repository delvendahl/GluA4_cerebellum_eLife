#pragma rtGlobals=3		// Use modern global access method and strict wave access.

// Simulation of GC with 4 Poisson MF inputs
// KO model
//
// Definition of model parameters:
//
// BEGIN definitions
static constant TIMESTEP = 0.05		// time step (ms)
static constant STIMSTART = 100		// start of stimulation (ms)
static constant DURATION = 1300		// simulation time (ms)
static constant SIMULATIONS = 10		// number of random trains per frequency
static constant DIAMETER = 11.7		// model diameter (µm), gives Cm=4.3005pF
static constant GTONIC = 0.0			// tonic GABA conductance (nS)
static constant GLEAK = 0.333			// leak conductance (nS)
static constant EGABA = -65 			// GABA reversal potential (mV)
static constant ELEAK = -110			// leak reversal potential (mV)
static constant VMEM = -80				// resting membrane potenital (mV)
static constant FREQ_INC = 20			// MF frequency increase (Hz)
static constant NUM_FREQ = 16			// number of simulations (running # SIMULATIONS each)
// END definitions

////////////////////////////////////
function init()

	NMTabAdd( "Pulse", "" )
	NMTabRemove( "Fit" )
	NMTabAdd( "Model", "" )
	NMTabRemove( "Event" )

	NMSet( tab="Model" )
	setrandomseed 0.4
	Set_model()

	NMSet( tab="Pulse" )
end

////////////////////////////////////
Function run()
	
	init()
	Make/O/D/N=(NUM_FREQ) MF_Freq, GC_Freq, GC_Freq_SD, Fano, VRE, Delay, Delay_SD
	MF_Freq = FREQ_INC + p * FREQ_INC

	variable i
	for(i = 0; i < NUM_FREQ; i += 1)
		Printf "Running model with: %.1f Hz MF stimulation\r", MF_Freq[i]

		wave results = run_model(MF_Freq[i])
		GC_Freq[i] = results[0]
		GC_Freq_SD[i] = results[1]
		Fano[i] = results[2]
		Delay[i] = results[3]
		Delay_SD[i] = results[4]

		string PrefStr
		sprintf PrefStr "P%04.0f_", MF_Freq[i]
		VRE[i] = GetVanRossumError(PrefStr)
		printf "\r\r"
	endfor

	KillWaves/Z results
	edit MF_Freq, GC_Freq, GC_Freq_SD, Fano, VRE, Delay, Delay_SD
	
	// Frequency
	Display GC_Freq vs MF_Freq
	execute "Graph1Style()"
	Label left "GC frequency (Hz)"

	// SD
	Display GC_Freq_SD vs MF_Freq
	execute "Graph1Style()"
	Label left "GC frequency SD (Hz)"

	// Fano factor
	Display Fano vs MF_Freq
	execute "Graph1Style()"
	Label left "Fano factor"

	// van Rossum error
	Display VRE vs MF_Freq
	execute "Graph1Style()"
	Label left "Van rossum distance"
	ModifyGraph log(left)=1

	// CV
	Duplicate/O GC_freq GC_freq_CV	
	GC_freq_CV = GC_Freq_SD/GC_Freq
	Display GC_Freq_CV vs MF_Freq
	execute "Graph1Style()"
	Label left "GC frequency CV"
	ModifyGraph log(left)=1
	
	// Delay
	Display Delay vs MF_Freq
	execute "Graph1Style()"
	Label left "First spike delay (ms)"
	
end

////////////////////////////////////
Function/WAVE run_model(frequency)
	variable frequency

	NMSet( tab="Pulse" )
	NM_make_Conduct_Waves(frequency)

	NMSet( tab="Model" )
	NMModelRun()

	NMSet( tab="Spike" )
	NMSpikeSet( xbgn=STIMSTART )
	NMSpikeSet( xend=DURATION - STIMSTART )
	NMSpikeRasterComputeAll( chanSelectList="A", waveSelectList="All", displayMode=1, delay=0 )
	NMSpikeRate( folder="Spike_Sim_Vmem_All_A", xRasterList="SP_RX_SimVmemAll_A0", xbgn=STIMSTART, xend=(DURATION - STIMSTART) )

	Wavestats/Q root:nmFolder0:Spike_Sim_Vmem_All_A:Rate_SP_RX_SimVmemAll_A0
	variable rate = V_avg * 1000
	variable rate_sd = V_sdev * 1000
	
	Wave SpikeDelays = CalcSpikeDelay()
	Wavestats/Q SpikeDelays
	variable delay = V_avg
	variable delay_sd = V_sdev
	Duplicate/O SpikeDelays, $("SpikeDelays_" +  num2str(frequency) + "Hz")
	
	Make/O/D/FREE/N=5 wOut = {rate, rate_sd, (rate_sd^2) / rate, delay, delay_sd}

	Printf "Simulation run with MF input at %.0f Hz\r",frequency
	Printf "GC frequency: %.5f ± %.5f Hz\r", rate, rate_sd
	Printf "First AP delay: %.5f ± %.5f ms\r", delay, delay_sd
	Printf "Fano factor: %.5f Hz\r",(rate_sd^2)/rate 

	NMSet( tab="Main" )

	string PrefStr
	sprintf PrefStr "P%04.0f_", frequency
	NMMainDuplicate( newPrefix=PrefStr )
	sort_into_folder(frequency)

	return wOut
End

////////////////////////////////////
Static Function Set_model()

	variable surface = DIAMETER^2*pi
	variable gTonicDensity = GTONIC / surface
	variable gLeakDensity = GLEAK / surface
	variable iClamp = GTONIC*(VMEM-EGABA) + GLEAK*(VMEM-ELEAK)
	
	NMSet( tab="Model" )	
	NMModelSet( modelSelect="IAF_AdEx" )
	NMModelVarSet( "SimulationTime", DURATION )
	NMModelVarSet( "TimeStep", TIMESTEP )
	NMModelVarSet( "Temperature", 37 )
	NMModelVarSet( "V0", VMEM )
	NMModelVarSet( "Diameter", DIAMETER )
	NMModelVarSet( "gTonicGABADensity", gTonicDensity )
	NMModelStrSet( "gGABA_WavePrefix", "gGABA_" )
	NMModelVarSet( "eGABA", EGABA )
	NMModelStrSet( "gNMDA_BlockFxnList", "GC_Schwartz2012" )
	NMModelVarSet( "gLeakDensity", gLeakDensity )
	NMModelVarSet( "eLeak", ELEAK )
	NMModelVarSet( "AP_Threshold", -45 )
	NMModelVarSet( "AP_Peak", 32 )
	NMModelVarSet( "AP_Reset", -75 )
	NMModelVarSet( "AP_Refrac", 0.75 )
	NMModelVarSet( "W_Tau", 75 )
	NMModelVarSet( "W_A", 0 )
	NMModelVarSet( "W_B", 0.35 )
	NMModelVarSet( "iClampAmp", iClamp )
	NMModelVarSet( "iClampOnset", 0 )
	NMModelVarSet( "iClampDuration", DURATION )
	NMModelVarSet( "LastSimulation", SIMULATIONS - 1 )
	NMModelTabUpdate()
	
End

////////////////////////////////////
Static Function sort_into_folder(frequency)
	Variable frequency

	DFREF saveDFR = GetDataFolderDFR()

	string PrefStr
	sprintf PrefStr "P%04.0f_", frequency

	NewDataFolder/O/S $PrefStr
	DFREF newDFR = GetDataFolderDFR()
	SetDataFolder saveDFR

	string wList = WaveList("Train*_"+num2istr(frequency) + "*", ";", "")
	wList += WaveList(PrefStr+"*", ";", "")
	wList += WaveList("SpikeDelays_*", ";", "")
	variable i
	for(i = 0; i < itemsinlist(wList); i += 1)
		wave w = $stringfromlist(i, wList)
		Movewave w, newDFR
	endfor
End

////////////////////////////////////
Function NM_make_Conduct_Waves(sumFreq) // code for generating MF-GC conductance waves
	Variable sumFreq // Hz, sum total freq of MF input

	sumFreq/= 1000 //kHz

	Variable error = -1
	Variable npnts, freq, interval, wcnt
	String tName0, tName1, tName2, tName3, wName, paramList, wList

	Variable SavePlasticityWavesAMPA = 0 // ( 0 ) no ( 1 ) yes
	Variable SavePlasticityWavesNMDA = 0

	Variable newRandomTrains = 1

	Variable numMFstim = 4 // simulating 4 MF inputs
	Variable refrac = 0.6 // ms

	Variable tbgn = STIMSTART
	Variable tend = DURATION - STIMSTART

	String gAMPA_Prefix = "gAMPA_"
	String gNMDA_Prefix = "gNMDA_"

	String fName = "nmFolder0"
	String df = "root:" + fName + ":"

	String directList1 = KO_NMPulseExp_GC_AMPAdirect( select = 0 ) + KO_NMPulseTrainRP_GC_AMPAdirect()
	String directList2 = KO_NMPulseExp_GC_AMPAdirect( select = 1 ) + KO_NMPulseTrainRP_GC_AMPAdirect()
	String spillList1 = KO_NMPulseExp_GC_AMPAspill( select = 0 ) + KO_NMPulseTrainRP_GC_AMPAspill()
	String spillList2 = KO_NMPulseExp_GC_AMPAspill( select = 1 ) + KO_NMPulseTrainRP_GC_AMPAspill()
	String spillList3 = KO_NMPulseExp_GC_AMPAspill( select = 2 ) + KO_NMPulseTrainRP_GC_AMPAspill()
	String nmdaList1 = WT_NMPulseExp_GC_NMDA( select = 0 ) + WT_NMPulseTrainRP_GC_NMDA()
	String nmdaList2 = WT_NMPulseExp_GC_NMDA( select = 1 ) + WT_NMPulseTrainRP_GC_NMDA()

	NMSet( folder="nmFolder0" )

	NMSet( tab="Pulse" )

	wList = WaveList( "PU_R_c*", ";", "" )
	NMKillWaves( wList )

	wList = WaveList( "PU_P_c*", ";", "" )
	NMKillWaves( wList )

	if ( newRandomTrains )
		wList = WaveList( "train*", ";", "" )
		NMKillWaves( wList )
	endif

	// create waves containing AP times

	freq = sumFreq / numMFstim
	interval = 1 / freq

	paramList = "train=random;interval=" + num2str( interval ) + ";tbgn=" + num2str( tbgn ) + ";tend=" + num2str( tend ) + ";refrac=" + num2str( refrac ) + ";"

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

	for ( wcnt = 0 ; wcnt < SIMULATIONS; wcnt += 1 )

		tName0 = "Train" + num2istr( wcnt ) + "_" + num2istr( sumFreq * 1000 ) + "Hz_0"
		tName1 = "Train" + num2istr( wcnt ) + "_" + num2istr( sumFreq * 1000 ) + "Hz_1"
		tName2 = "Train" + num2istr( wcnt ) + "_" + num2istr( sumFreq * 1000 ) + "Hz_2"
		tName3 = "Train" + num2istr( wcnt ) + "_" + num2istr( sumFreq * 1000 ) + "Hz_3"

		if ( newRandomTrains || !WaveExists( $tName0 ) )
			if ( NMPulseTrainRandomTimes( df, tName0, paramList, DURATION ) <= 0 )
				//return error
			endif
			if ( NMPulseTrainRandomTimes( df, tName1, paramList, DURATION ) <= 0 )
				//return error
			endif
			if ( NMPulseTrainRandomTimes( df, tName2, paramList, DURATION ) <= 0 )
				//return error
			endif
			if ( NMPulseTrainRandomTimes( df, tName3, paramList, DURATION ) <= 0 )
				//return error
			endif
		endif

		// gAMPA (direct + spillover)
		NMPulseSet( wavePrefix=gAMPA_Prefix, update=1 )
		NMPulseSet( wavePrefix=gAMPA_Prefix, update=1 )

		NMPulseConfigAdd( "wave=" + num2istr( wcnt ) + ";train=" + tName0 + ";" + directList1 )
		NMPulseConfigAdd( "wave=" + num2istr( wcnt ) + ";train=" + tName0 + ";" + directList2 )
		NMPulseConfigAdd( "wave=" + num2istr( wcnt ) + ";train=" + tName0 + ";" + spillList1 )
		NMPulseConfigAdd( "wave=" + num2istr( wcnt ) + ";train=" + tName0 + ";" + spillList2 )
		NMPulseConfigAdd( "wave=" + num2istr( wcnt ) + ";train=" + tName0 + ";" + spillList3 )

	  	NMPulseConfigAdd( "wave=" + num2istr( wcnt ) + ";train=" + tName1 + ";" + directList1 )
	  	NMPulseConfigAdd( "wave=" + num2istr( wcnt ) + ";train=" + tName1 + ";" + directList2 )
	  	NMPulseConfigAdd( "wave=" + num2istr( wcnt ) + ";train=" + tName1 + ";" + spillList1 )
	  	NMPulseConfigAdd( "wave=" + num2istr( wcnt ) + ";train=" + tName1 + ";" + spillList2 )
	  	NMPulseConfigAdd( "wave=" + num2istr( wcnt ) + ";train=" + tName1 + ";" + spillList3 )

	  	NMPulseConfigAdd( "wave=" + num2istr( wcnt ) + ";train=" + tName2 + ";" + directList1 )
	  	NMPulseConfigAdd( "wave=" + num2istr( wcnt ) + ";train=" + tName2 + ";" + directList2 )
	  	NMPulseConfigAdd( "wave=" + num2istr( wcnt ) + ";train=" + tName2 + ";" + spillList1 )
	  	NMPulseConfigAdd( "wave=" + num2istr( wcnt ) + ";train=" + tName2 + ";" + spillList2 )
	  	NMPulseConfigAdd( "wave=" + num2istr( wcnt ) + ";train=" + tName2 + ";" + spillList3 )

	  	NMPulseConfigAdd( "wave=" + num2istr( wcnt ) + ";train=" + tName3 + ";" + directList1 )
	  	NMPulseConfigAdd( "wave=" + num2istr( wcnt ) + ";train=" + tName3 + ";" + directList2 )
	  	NMPulseConfigAdd( "wave=" + num2istr( wcnt ) + ";train=" + tName3 + ";" + spillList1 )
	  	NMPulseConfigAdd( "wave=" + num2istr( wcnt ) + ";train=" + tName3 + ";" + spillList2 )
	  	NMPulseConfigAdd( "wave=" + num2istr( wcnt ) + ";train=" + tName3 + ";" + spillList3 )

		// gNMDA
		NMPulseSet(  wavePrefix=gNMDA_Prefix )
		NMPulseSet(  wavePrefix=gNMDA_Prefix )

		NMPulseConfigAdd( "wave=" + num2istr( wcnt ) + ";train=" + tName0 + ";" + nmdaList1)
		NMPulseConfigAdd( "wave=" + num2istr( wcnt ) + ";train=" + tName0 + ";" + nmdaList2 )
	  	NMPulseConfigAdd( "wave=" + num2istr( wcnt ) + ";train=" + tName1 + ";" + nmdaList1 )
	  	NMPulseConfigAdd( "wave=" + num2istr( wcnt ) + ";train=" + tName1 + ";" + nmdaList2 )
	  	NMPulseConfigAdd( "wave=" + num2istr( wcnt ) + ";train=" + tName2 + ";" + nmdaList1 )
	  	NMPulseConfigAdd( "wave=" + num2istr( wcnt ) + ";train=" + tName2 + ";" + nmdaList2 )
	  	NMPulseConfigAdd( "wave=" + num2istr( wcnt ) + ";train=" + tName3 + ";" + nmdaList1 )
	  	NMPulseConfigAdd( "wave=" + num2istr( wcnt ) + ";train=" + tName3 + ";" + nmdaList2 )

	endfor

	NMPulseSet( wavePrefix=gAMPA_Prefix, update=1 )
	NMPulseSet( wavePrefix=gAMPA_Prefix, update=1 )
	NMConfigVarSet( "Pulse", "SavePlasticityWaves", SavePlasticityWavesAMPA )
	NMPulseExecute()


	NMPulseSet( wavePrefix=gNMDA_Prefix, update=1 )
	NMPulseSet( wavePrefix=gNMDA_Prefix, update=1 )
	NMConfigVarSet( "Pulse", "SavePlasticityWaves", SavePlasticityWavesNMDA )
	NMPulseExecute()

End

////////////////////////////////////
// functions generating parameter lists for AMPA & NMDA
////////////////////////////////////
// CONDUCTANCE FUNCTIONS FOR WT
////////////////////////////////////
Static Function /S WT_NMPulseExp_GC_AMPAdirect( [ select, p ] ) // Billings 2014
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
			amp = 1.4 // nS
			pp.tau2 =  0.5 // ms
			break
		case 1:
			amp = 0.55 // nS
			pp.tau2 = 4 // ms
			break
	endswitch

	pp.amp3 = 0
	pp.amp4 = 0

	return "pulse=exp;amp=" + num2str( amp ) + ";" + NMPulseExpParamList( pp )

End // WT_NMPulseExp_GC_AMPAdirect

////////////////////////////////////
Static Function /S WT_NMPulseTrainRP_GC_AMPAdirect( [ p ] )
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

////////////////////////////////////
Static Function /S WT_NMPulseExp_GC_AMPAspill( [ select, p ] ) // Billings 2014
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
			amp = 0.1 // nS
			pp.tau2 = 0.8 // ms
			break
		case 1:
			amp = 0.25 // nS
			pp.tau2 = 4 // ms
			break
		case 2:
			amp = 0.15 // nS
			pp.tau2 = 30 // ms
			break
		default:
			return ""
	endswitch

	return "pulse=exp;amp=" + num2str( amp ) + ";" + NMPulseExpParamList( pp )

End // WT_NMPulseExp_GC_AMPAspill

////////////////////////////////////
Static Function /S WT_NMPulseTrainRP_GC_AMPAspill( [ p ] )
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

////////////////////////////////////
Static Function /S WT_NMPulseExp_GC_NMDA( [ select, p ] ) // Billings 2014
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
			amp = 0.3 // nS
			pp.tau2 = 30 // ms
			break
		case 1:
			amp = 0.14 // nS
			pp.tau2 = 70 // ms
			break
		default:
			return ""
	endswitch

	return "pulse=exp;amp=" + num2str( amp ) + ";" + NMPulseExpParamList( pp )

End // WT_NMPulseExp_GC_NMDA

////////////////////////////////////
Static Function /S WT_NMPulseTrainRP_GC_NMDA( [ p ] )
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

//////////////////////////////////
// CONDUCTANCE FUNCTIONS FOR KO
//////////////////////////////////
Static Function /S KO_NMPulseExp_GC_AMPAdirect( [ select, p ] ) // Billings 2014
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
			amp = 0.196 // nS
			pp.tau2 =  0.9 // ms
			break
		case 1:
			amp = 0.077 // nS
			pp.tau2 = 7.2 // ms
			break
	endswitch

	return "pulse=exp;amp=" + num2str( amp ) + ";" + NMPulseExpParamList( pp )
End // KO_NMPulseExp_GC_AMPAdirect

////////////////////////////////////
Static Function /S KO_NMPulseTrainRP_GC_AMPAdirect( [ p ] )
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

////////////////////////////////////
Static Function /S KO_NMPulseExp_GC_AMPAspill( [ select, p ] ) // Billings 2014
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
			amp = 0.014 // nS
			pp.tau2 = 1.44 // ms
			break
		case 1:
			amp = 0.035 // nS
			pp.tau2 = 7.2 // ms
			break
		case 2:
			amp = 0.021 // nS
			pp.tau2 = 54 // ms
			break
		default:
			return ""
	endswitch

	return "pulse=exp;amp=" + num2str( amp ) + ";" + NMPulseExpParamList( pp )
End // KO_NMPulseExp_GC_AMPAspill

////////////////////////////////////
Static Function /S KO_NMPulseTrainRP_GC_AMPAspill( [ p ] )
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

////////////////////////////////////
Static Function /S KO_NMPulseExp_GC_NMDA( [ select, p ] ) // Billings 2014
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
			amp = 0.45 // nS
			pp.tau2 = 30 // ms
			break
		case 1:
			amp = 0.21 // nS
			pp.tau2 = 70 // ms
			break
		default:
			return ""
	endswitch

	return "pulse=exp;amp=" + num2str( amp ) + ";" + NMPulseExpParamList( pp )
End // KO_NMPulseExp_GC_NMDA

////////////////////////////////////
Static Function /S KO_NMPulseTrainRP_GC_NMDA( [ p ] )
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

////////////////////////////////////
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

////////////////////////////////////
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

////////////////////////////////////
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

////////////////////////////////////
Static Function/WAVE CalcSpikeDelay()
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

////////////////////////////////////
Proc Graph1Style() : GraphStyle
	PauseUpdate; Silent 1		// modifying window...
	ModifyGraph/Z expand=2,gFont="Verdana",gfSize=8,width=85.0394,height=85.0394
	ModifyGraph/Z mode=4
	ModifyGraph/Z marker=19
	ModifyGraph/Z rgb=(1,34817,52428)
	ModifyGraph/Z nticks=3
	ModifyGraph/Z lblMargin=5
	ModifyGraph/Z standoff=0
	ModifyGraph/Z axThick=0.75
	ModifyGraph/Z btLen=3
	ModifyGraph/Z tlOffset(bottom)=-2
	Label/Z left "GC frequency SD (Hz)"
	Label/Z bottom "Sum MF frequency (Hz)"
	SetAxis/Z/A/N=1/E=1 left
	SetAxis/Z/A/E=1 bottom
EndMacro