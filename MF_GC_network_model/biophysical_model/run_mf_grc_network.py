# Runs biophysical simulations of MF-GC network
# Designed to write files to be run in parallel on a server or cluster, 
# but can also be used to run individually
# python run_mf_grc_network.py
# Generates xml and nml files for model, saved in tempdata folder
# To simulate, switch run to True, spike times are saved as .dat

# Flow for biophysical model:
# python initialize_network_params.py to generate params_file.pkl for a specific correlation radius (sigma)
# python run_mf_grc_network.py
# python save_samples_as_txt.py basedir to convert .dat files of spiketimes to .txt files of activity patterns
# Then can run run_learning.py, get_spar_cov.py, etc.

import neuroml as nml
from pyneuroml import pynml
from pyneuroml.lems.LEMSSimulation import LEMSSimulation
import lems.api as lems
import numpy as np
import os
import pickle as pkl
import scipy.io as io
import random
from datetime import datetime
from multiprocessing import Condition, Pool
from functools import partial


def generate_grc_layer_network(runID, correlationRadius, duration=180.0, dt=0.05,
                               minimumISI=2.0, ONRate=50.0, OFFRate=0.0, condition='orig', run=False):
    """
    Creates GrC layer network and runs LEMS simulation.

    :param runID: ID number of the run (=pattern)
    :param correlationRadius: MF correlation radius, out of [0,5,10,15,20,25,30]
    :param duration: total duration of the simulation (default=180ms)
    :param dt: simulation time step (default=0.05ms)
    :param minimumISI: minimum ISI for MF (ms, default=2.0)
    :param ONRate: rate of active MFs (Hz, default=50.0)
    :param OFFRate: rate of inactive MFs (Hz, default=0.0)
    :param condition: which simulation to run ('orig' or 'ko', default='orig')
    :param run: boolean (default=False)
    :return: boolean (success)
    """

    # Load parameters for this run from file
    file = open('../params_file.pkl', 'r')
    p = pkl.load(file)
    N_syn = p['N_syn'][int(runID)]
    f_mf = p['f_mf'][int(runID)]
    run_num = p['run_num'][int(runID)]
    file.close()

    # set the random number generator seed
    np.random.seed(runID)
    random.seed(runID)
 
    # get connectivity matrix between cells from file
    file = open('../../network_structures/GCLconnectivity_{:.0f}.pkl'.format(N_syn))
    p = pkl.load(file)
    file.close()
    conn_mat = p['conn_mat']
    N_mf, N_grc = conn_mat.shape
    assert (np.all(conn_mat.sum(axis=0) == N_syn)), 'Connectivity matrix is incorrect.'

    # turn fraction of MFs off
    if correlationRadius == 0:  # Activate MFs randomly
        N_mf_ON = int(N_mf * f_mf)
        mf_indices_ON = random.sample(range(N_mf), N_mf_ON)
        mf_indices_ON.sort()
    elif correlationRadius > 0:  # Spatially correlated MFs
        f_mf_range = np.linspace(.05, .95, 19)
        f_mf_ix = np.where(np.isclose(f_mf_range, f_mf))[0][0]
        p = io.loadmat('../../input_statistics/mf_patterns_r{:.0f}.mat'.format(correlationRadius))
        R = p['Rs'][:, :, f_mf_ix]
        g = p['gs'][f_mf_ix]
        t = np.dot(R.transpose(), np.random.randn(N_mf))
        S = (t > -g * np.ones(N_mf))
        mf_indices_ON = np.where(S)[0]
        N_mf_ON = len(mf_indices_ON)

    N_mf_OFF = N_mf - N_mf_ON
    mf_indices_OFF = [x for x in range(N_mf) if x not in mf_indices_ON]
    mf_indices_OFF.sort()

    # load NeuroML components, LEMS components and LEMS componentTypes from external files
    # Spike generator (for Poisson MF spiking)
    spike_generator_file_name = "../../grc_lemsDefinitions/spikeGenerators.xml"
    spike_generator_doc = pynml.read_lems_file(spike_generator_file_name)

    # Integrate-and-fire GC model
    iaf_nml2_file_name = '../../grc_lemsDefinitions/IaF_GrC_{}.nml'.format(condition)
    iaF_GrC_doc = pynml.read_neuroml2_file(iaf_nml2_file_name)
    iaF_GrC = iaF_GrC_doc.iaf_ref_cells[0]

    # AMPAR and NMDAR mediated synapses
    ampa_syn_filename = "../../grc_lemsDefinitions/MFGrC_AMPA_{}.xml".format(condition)
    nmda_syn_filename = "../../grc_lemsDefinitions/MFGrC_NMDA_{}.xml".format(condition)
    rothmanMFToGrCAMPA_doc = pynml.read_lems_file(ampa_syn_filename)
    rothmanMFToGrCNMDA_doc = pynml.read_lems_file(nmda_syn_filename)
    
    # Define components from the componentTypes we just loaded
    # Refractory poisson input -- representing active MF
    spike_generator_ref_poisson_type = spike_generator_doc.component_types['MyspikeGeneratorRefPoisson']
    lems_instances_doc = lems.Model()
    spike_generator_on = lems.Component("mossySpikerON", spike_generator_ref_poisson_type.name)
    spike_generator_on.set_parameter("minimumISI", "%s ms" % minimumISI)
    spike_generator_on.set_parameter("averageRate", "%s Hz" % ONRate)
    lems_instances_doc.add(spike_generator_on)

    # Refractory poisson input -- representing silent MF
    spike_generator_off = lems.Component("mossySpikerOFF", spike_generator_ref_poisson_type.name)
    spike_generator_off.set_parameter("minimumISI", "%s ms" % minimumISI)
    spike_generator_off.set_parameter("averageRate", "%s Hz" % OFFRate)
    lems_instances_doc.add(spike_generator_off)

    # Synapses
    rothmanMFToGrCAMPA = rothmanMFToGrCAMPA_doc.components['RothmanMFToGrCAMPA'].id
    rothmanMFToGrCNMDA = rothmanMFToGrCNMDA_doc.components['RothmanMFToGrCNMDA'].id
    
    # Create ON MF, OFF MF, and GC populations
    GrCPop = nml.Population(id="GrCPop", component=iaF_GrC.id, size=N_grc)
    mossySpikersPopON = nml.Population(id=spike_generator_on.id + "Pop", component=spike_generator_on.id, size=N_mf_ON)
    mossySpikersPopOFF = nml.Population(id=spike_generator_off.id + "Pop", component=spike_generator_off.id,
                                        size=N_mf_OFF)
    
    # Create network and add populations
    net = nml.Network(id="network")
    net_doc = nml.NeuroMLDocument(id=net.id)
    net_doc.networks.append(net)
    net.populations.append(GrCPop)
    net.populations.append(mossySpikersPopON)
    net.populations.append(mossySpikersPopOFF)
    
    # MF-GC connectivity
    # First connect ON MFs to GCs
    for mf_ix_ON in range(N_mf_ON):
        mf_ix = mf_indices_ON[mf_ix_ON]
        # Find which GCs are neighbors
        innervated_grcs = np.where(conn_mat[mf_ix, :] == 1)[0]
        for grc_ix in innervated_grcs:
            # Add AMPAR and NMDAR mediated synapses
            for synapse in [rothmanMFToGrCAMPA, rothmanMFToGrCNMDA]:
                connection = nml.SynapticConnection(from_='{}[{}]'.format(mossySpikersPopON.id, mf_ix_ON),
                                                    synapse=synapse,
                                                    to='GrCPop[{}]'.format(grc_ix))
                net.synaptic_connections.append(connection)
    
    # Now connect OFF MFs to GCs
    for mf_ix_OFF in range(N_mf_OFF):
        mf_ix = mf_indices_OFF[mf_ix_OFF]
        # Find which GCs are neighbors
        innervated_grcs = np.where(conn_mat[mf_ix, :] == 1)[0]
        for grc_ix in innervated_grcs:
            # Add AMPAR and NMDAR mediated synapses
            for synapse in [rothmanMFToGrCAMPA, rothmanMFToGrCNMDA]:
                connection = nml.SynapticConnection(from_='{}[{}]'.format(mossySpikersPopOFF.id, mf_ix_OFF),
                                                    synapse=synapse,
                                                    to='GrCPop[{}]'.format(grc_ix))
                net.synaptic_connections.append(connection)
    
    # Write network to file
    net_file_name = 'generated_network_{}.net.nml'.format(runID)
    pynml.write_neuroml2_file(net_doc, net_file_name, validate=False)

    # Write LEMS instances to file
    lems_instances_file_name = 'instances_{}.xml'.format(runID)
    pynml.write_lems_file(lems_instances_doc, lems_instances_file_name, validate=False)

    # Create a LEMSSimulation to manage creation of LEMS file
    ls = LEMSSimulation('sim_{}'.format(runID), duration, dt)

    # Point to network as target of simulation
    ls.assign_simulation_target(net.id)

    # Include generated/existing NeuroML2 files
    ls.include_neuroml2_file(iaf_nml2_file_name)
    ls.include_lems_file(spike_generator_file_name, include_included=False)
    ls.include_lems_file(lems_instances_file_name)
    ls.include_lems_file(ampa_syn_filename, include_included=False)
    ls.include_lems_file(nmda_syn_filename, include_included=False)
    ls.include_neuroml2_file(net_file_name)

    # Specify Displays and Output Files
    # Details for saving output files
    basedir = '../results/{}_data_r{}/'.format(condition, correlationRadius)

    # Add parameter values to spike time filename
    end_filename = '{}_{:.2f}_{}'.format(N_syn, f_mf, run_num)

    # Save MF spike times under basedir + MF_spikes_ + end_filename
    eof0 = 'MFspikes_file'
    ls.create_event_output_file(eof0, basedir + "MF_spikes_" + end_filename + ".dat")
    # ON MFs
    for i in range(mossySpikersPopON.size):
        ls.add_selection_to_event_output_file(eof0, mf_indices_ON[i], "%s[%i]" % (mossySpikersPopON.id, i), 'spike')
    # OFF MFs
    for i in range(mossySpikersPopOFF.size):
        ls.add_selection_to_event_output_file(eof0, mf_indices_OFF[i], "%s[%i]" % (mossySpikersPopOFF.id, i), 'spike')
    
    # Save GC spike times under basedir + GrC_spikes_ + end_filename
    eof1 = 'GrCspikes_file'
    ls.create_event_output_file(eof1, basedir + "GrC_spikes_" + end_filename + ".dat")
    for i in range(GrCPop.size):
        ls.add_selection_to_event_output_file(eof1, i, "%s[%i]" % (GrCPop.id, i), 'spike')
    
    lems_file_name = ls.save_to_file()
    
    if run:
        results = pynml.run_lems_with_jneuroml(lems_file_name,
                                               max_memory="8G",
                                               nogui=True,
                                               load_saved_data=False,
                                               plot=False)

        return results


if __name__ == '__main__':
    
    startTime = datetime.now()
    
    # type of simulation to run; change this to 'ko' for GluA4_KO model
    sim_type = 'orig' 

    # total number of simulation runs per correlation radius (= num_patterns * len(f_mf))
    runs = range(5760)

    # run simulation for the following MF correlation radii
    corrs = [0,5,10,15,20,25,30]

    # Move working dir to tempdata to hide xml and nml files
    os.chdir('tempdata')
    
    pool = Pool()
    for i, c in enumerate(corrs):
        run_this = partial(generate_grc_layer_network,
                           correlationRadius=c,
                           duration=180,
                           dt=0.05,
                           minimumISI=2,
                           ONRate=50,
                           OFFRate=0,
                           condition=sim_type,
                           run=True)
        results = pool.map(run_this, runs)
    
    pool.close()
    pool.join()
    
    print(datetime.now() - startTime)
    
    # Move back to original directory
    os.chdir('..')
