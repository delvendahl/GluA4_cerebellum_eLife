#pragma rtGlobals=3		// Use modern global access method and strict wave access.

// Simulation of cerebellar GC
// Runs integrate-and-fire model with 2pA current step increase
// Simulates control with ("ctrl") and without tonic inhibition ("bmi")


// Definition of model parameters:
//
// BEGIN definitions
static constant TIMESTEP = 0.05	// time step (ms)
static constant STIMSTART = 50		// start of stimulation (ms)
static constant STIMDURATION = 200	// duration of current injection
static constant DURATION = 300		// simulation time (ms)
static constant SIMULATIONS = 31	// number of simulations
static constant CURRENT_STEP = 2	// current increase step (pA)
static constant DIAMETER = 11.7	// model diameter (um), gives Cm=4.3005pF
static constant RMP = -100			// resting membrane potential (mV)
static constant TEMP = 25			// temperature (deg C)
static constant GTONIC = 0.16		// tonic GABA conductance (nS)
static constant GLEAK = 0.366		// leak conductance (nS)
static constant EGABA = -65 		// GABA reversal potential (mV)
static constant ELEAK = -110		// leak reversal potential (mV)
// END definitions


// run this function
function run()

	init()
	
	run_ctrl()
	run_bmi()
	
end


// initialize Nueromatic and set random seed
static function init()
	NMTabAdd( "Pulse", "" )
	NMTabRemove( "Fit" )
	NMTabAdd( "Model", "" )
	NMTabRemove( "Event" )
	setrandomseed 0.4
end


// run control ("ctrl")
static function run_ctrl()
		
	variable surface = DIAMETER^2*pi
	variable gTonicDensity = GTONIC / surface
	variable gLeakDensity = GLEAK / surface
	variable iClamp = GTONIC*(RMP-EGABA) + GLEAK*(RMP-ELEAK)
	
	NMSet( tab="Pulse" )
	NMPulseSet( dx=TIMESTEP )
	NMPulseSet( waveLength=DURATION )
	
	String paramList = "wave=all;pulse=square;" + "amp=" + num2str( - iClamp / RMP ) + ";onset=0" + ";width=" + num2str( DURATION ) + ";"
	
	NMPulseConfigRemove( all=1 )
	NMPulseConfigAdd( paramList )
	NMPulseSet( numWaves=SIMULATIONS )
	NMPulseExecute()
	
	NMSet( tab="Model" )
	Set_Model()
	NMModelStrSet( "WavePrefix", "ctrl_" )
	NMModelVarSet( "eLeak", ELEAK )
	
	NMModelStrSet( "gAMPA_WavePrefix", "Pulse" )
	NMModelVarSet( "eAMPA", 0 )
	
	NMModelVarSet( "iClampAmp", 0 )
	NMModelVarSet( "iClampAmpInc", CURRENT_STEP )
	
	NMModelVarSet( "gTonicGABADensity", gTonicDensity )
	NMModelVarSet( "gLeakDensity", gLeakDensity )
	NMModelVarSet( "eLeak", ELEAK )

	NMModelTabUpdate()
	
	NMModelRun()
	
	NMSet( tab="Spike" )
	NMSpikeSet( xbgn=STIMSTART )
	NMSpikeSet( xend=STIMSTART+STIMDURATION )
	NMSpikeRasterComputeAll( chanSelectList="A", waveSelectList="All", displayMode=1, delay=0 )
	
	NMSpikeRate( folder="Spike_ctrl_Vmem_All_A", xRasterList="SP_RX_ctrlVmemAll_A0", xbgn=STIMSTART, xend=STIMSTART+STIMDURATION )
	
end


// run control without inhibition ("bmi")
static function run_bmi()
	
	variable surface = DIAMETER^2*pi
	variable gTonicDensity = 0.0		// no tonic inhibition
	variable gLeakDensity = GLEAK / surface
	variable iClamp =  GLEAK*(RMP-ELEAK) // no tonic inhibition
	NMSet( tab="Pulse" )
	NMPulseSet( dx=TIMESTEP )
	NMPulseSet( waveLength=DURATION )
	
	String paramList = "wave=all;pulse=square;" + "amp=" + num2str( - iClamp /RMP ) + ";onset=0" + ";width=" + num2str( DURATION ) + ";"
	
	NMPulseConfigRemove( all=1 )
	NMPulseConfigAdd( paramList )
	NMPulseSet( numWaves=SIMULATIONS )
	NMPulseExecute()

	NMSet( tab="Model" )
	Set_Model()
	NMModelStrSet( "WavePrefix", "bmi_" )
	NMModelVarSet( "eLeak", ELEAK )
	
	NMModelStrSet( "gAMPA_WavePrefix", "Pulse" )
	NMModelVarSet( "eAMPA", 0 )
	
	NMModelVarSet( "iClampAmp", 0 )
	NMModelVarSet( "iClampAmpInc", CURRENT_STEP )
	
	NMModelVarSet( "gTonicGABADensity", gTonicDensity )
	NMModelVarSet( "gLeakDensity", gLeakDensity )
	NMModelVarSet( "eLeak", ELEAK )

	NMModelTabUpdate()
	
	NMModelRun()
	
	NMSet( tab="Spike" )
	NMSpikeSet( xbgn=STIMSTART )
	NMSpikeSet( xend=STIMSTART+STIMDURATION )
	NMSpikeRasterComputeAll( chanSelectList="A", waveSelectList="All", displayMode=1, delay=0 )
	
	NMSpikeRate( folder="Spike_bmi_Vmem_All_A", xRasterList="SP_RX_bmiVmemAll_A0", xbgn=STIMSTART, xend=STIMSTART+STIMDURATION )
	
end


// set model parameters
static function Set_Model()
	NMSet( tab="Model" )

	NMModelSet( modelSelect="IAF_AdEx" )
	
	NMModelVarSet( "Diameter", DIAMETER )
	NMModelVarSet( "LastSimulation", SIMULATIONS-1 )
	NMModelVarSet( "SimulationTime", DURATION )
	NMModelVarSet( "TimeStep", TIMESTEP )
	
	NMModelVarSet( "Temperature", TEMP )
	NMModelVarSet( "V0", RMP )	
	NMModelVarSet( "iClampAmp", 0 )
	NMModelVarSet( "iClampAmpInc", 2 )
	NMModelVarSet( "iClampOnset", STIMSTART )
	NMModelVarSet( "iClampDuration", STIMDURATION )
	
	NMModelVarSet( "eGABA", EGABA )
	
	NMModelVarSet( "AP_Threshold", -45 )
	NMModelVarSet( "AP_ThreshSlope", 2 )
	NMModelVarSet( "AP_Peak", 32 )
	NMModelVarSet( "AP_Reset", -75 )
	NMModelVarSet( "AP_Refrac", 1.2 )
	NMModelVarSet( "W_Tau", 75 )
	NMModelVarSet( "W_A", 0 )
	NMModelVarSet( "W_B", 0.25 )
	
	NMModelTabUpdate()
end
