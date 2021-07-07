import matplotlib
matplotlib.use('TkAgg')
from pyneuroml.analysis import generate_current_vs_frequency_curve
import os


if __name__ == '__main__':

    sim_type = 'orig' # change to 'ko' for GluA4_KO model
    
    os.chdir('tempdata')
    cellname = '../../grc_lemsDefinitions/IaF_GrC_{}.nml'.format(sim_type)

    generate_current_vs_frequency_curve(cellname,
                                        'IaF_GrC',
                                        start_amp_nA=0.0,
                                        end_amp_nA=0.1,
                                        step_nA=0.005,
                                        analysis_duration=200,
                                        analysis_delay=0,
                                        pre_zero_pulse=50,
                                        post_zero_pulse=50,
                                        spike_threshold_mV=-46,
                                        plot_voltage_traces=False,
                                        plot_if=True,
                                        plot_iv=False,
                                        save_if_data_to = '../results/IF_result_{}.txt'.format(sim_type))

    # Move back to original directory
    os.chdir('..')
