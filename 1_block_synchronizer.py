import pypulseq as puls
import scipy as sp
import numpy as np
import json
from matplotlib import pyplot as plt

data_seq = sp.io.loadmat('Sequence_from_seq_SE_rfdeath_5000.mat')
A_data = data_seq['A']
seq_file = "SE_rfdeath_5000.seq"
seq_input = puls.Sequence()
seq_input.read(file_path=seq_file)
seq_output_dict = seq_input.waveforms_export(time_range=(0, 3))

# added type check in Sequence.block, read does not make an empty variable with a type
# is there the other way to do it?
# print(seq_output_dict['gx'])
# Engage what exactly every array means
# print(seq_input.waveforms_and_times())
# plt.plot()
# plt.show()

# print(seq_output_dict)
# t_adc t_rf t_rf_centers t_gx t_gy t_gz adc rf rf_centers gx gy gz
# seq_input.plot()

# plt.plot(seq_output_dict['t_rf'], seq_output_dict['rf'])
# plt.show()
# plt.plot(seq_output_dict['t_adc'], seq_output_dict['adc'])
# plt.show()


local_definitions = seq_input.definitions
ADC_raster = local_definitions['AdcRasterTime']
RF_raster = local_definitions['RadiofrequencyRasterTime']

RF_dtime = 1000 * 1e-6
TR_dtime = 1000 * 1e-6
# artificial delays

time_info = seq_input.duration()
blocks_number = time_info[1]
time_dur = time_info[0]
time_step = 1e-6
round_const_rf = 7
round_const_adc = 5
# TODO: how to auto round_const
N_samples = int(time_dur / time_step)
# TODO: why two times bigger? what effort on output
time_sample = np.linspace(0, time_dur, N_samples)

gate_adc = np.zeros(N_samples)
gate_rf = np.zeros(N_samples)
gate_tr_switch = np.ones(N_samples)
gate_gx = np.zeros(N_samples)
gate_gy = np.zeros(N_samples)
gate_gz = np.zeros(N_samples)

local_delay_rf = RF_dtime
local_delay_tr = TR_dtime
local_raster_time = time_step

# TODO: function defining beginning and ending of the RF events
RF_assintant = [seq_output_dict['t_rf'][0] - RF_dtime, seq_output_dict['t_rf'][-1]]


def gates_output(gates):
    with open('data_output/tr_switch.txt', 'w') as f:
        data = str(tuple(gates['tr_switch']))
        f.write(data)
    with open('data_output/rf.txt', 'w') as f:
        data = str(tuple(gates['rf']))
        f.write(data)
    with open('data_output/adc.txt', 'w') as f:
        data = str(tuple(gates['adc']))
        f.write(data)
    data = {'gate_gx': tuple(gates['gx']),
            'gate_gy': tuple(gates['gy']),
            'gate_gz': tuple(gates['gz'])}
    with open('data_output/gradient_gates.json', 'w') as outfile:
        json.dump(data, outfile)


def adc_correction():
    rise_time, fall_time = None, None
    is_adc_inside = False
    for j in range(blocks_number - 1):
        iterable_block = seq_input.get_block(block_index=j + 1)
        if iterable_block.adc is not None:
            is_adc_inside = True
            rise_time = iterable_block.gx.rise_time
            fall_time = iterable_block.gx.fall_time
    if not is_adc_inside:
        raise Exception("No ADC event found inside sequence")
    return rise_time, fall_time


def adc_event_edges(local_gate_adc):
    num_begin_l = 0
    flag_begin = False
    flag_finish = False
    num_finish_l = 1
    for k in range(len(local_gate_adc) - 1):
        if local_gate_adc[k] != 0 and not flag_begin:
            num_begin_l = k
            flag_begin = True
        if local_gate_adc[k] != 0 and local_gate_adc[k + 1] == 0 and not flag_finish:
            num_finish_l = k
            flag_finish = True
    return num_begin_l, num_finish_l


for i in range(N_samples):
    # delaying of RF event for time period of local delay
    if RF_assintant[0] - RF_raster < time_sample[i] < RF_assintant[0] + RF_raster:
        RF_stop = int(RF_assintant[1] / time_step)
        gate_rf[i:RF_stop] = 1.0
        var = 1
    # mandatory disabling of RF gate due to ADC work same time
    gate_rf_2 = map(lambda x: time_sample[i] - ADC_raster < x < time_sample[i] + ADC_raster and 1 or 0,
                    seq_output_dict['t_adc'])
    if np.any(np.array(list(gate_rf_2)) > 0):
        gate_rf[i] = 0.0
    # TR switch with own delay before ADC turning
    gate_tr_1 = map(lambda x: time_sample[i] - ADC_raster < x < time_sample[i] + ADC_raster and 1 or 0,
                    seq_output_dict['t_adc'])
    if np.any(np.array(list(gate_tr_1)) > 0):
        block_delay_tr = int(local_delay_tr / time_step)
        gate_tr_switch[i - block_delay_tr:i + 1] = 0.0
    # first step of ADC gate - enabling
    # TODO: ADC gate feeling gradients form of rise and fall
    gate_adc_1 = map(lambda x: time_sample[i] - ADC_raster < x < time_sample[i] + ADC_raster and 1 or 0,
                     seq_output_dict['t_adc'])
    if np.any(np.array(list(gate_adc_1)) > 0):
        gate_adc[i] = 1.0

# adc correction sue to rise and fall time of gradient
# defining time that ADC need to be disabled during of
rise_time_loc, fall_time_loc = adc_correction()
num_beg, num_fin = adc_event_edges(gate_adc)
rise_time_tick = int(rise_time_loc / time_step)
fall_time_tick = int(rise_time_loc / time_step)
gate_adc[num_beg:num_beg + rise_time_tick] = 0.0
gate_adc[num_fin - fall_time_tick:num_fin+1] = 0.0

gates_release = {"adc": gate_adc,
                 "rf": gate_rf,
                 "tr_switch": gate_tr_switch,
                 "gx": gate_gx,
                 "gy": gate_gy,
                 "gz": gate_gz}

plt.plot(seq_output_dict['t_gx'][:int(N_samples)], seq_output_dict['gx'][:int(N_samples)])
plt.plot(seq_output_dict['t_gy'][:int(N_samples)], seq_output_dict['gy'][:int(N_samples)])
plt.plot(seq_output_dict['t_gz'][:int(N_samples)], seq_output_dict['gz'][:int(N_samples)])
plt.show()

plt.plot(seq_output_dict['t_gx'][:int(N_samples)], seq_output_dict['gx'][:int(N_samples)] / 720)
plt.plot(time_sample[:int(N_samples / 3)], gate_adc[:int(N_samples / 3)], label='ADC gate')
# plt.plot(time_sample[:int(N_samples / 3)], gate_tr_switch[:int(N_samples / 3)], label='TR switch')
# plt.plot(seq_output_dict['t_rf'], seq_output_dict['rf'] / 210, label='RF signal')
# plt.plot(time_sample[:int(N_samples / 3)], gate_rf[:int(N_samples / 3)], label='RF gate')
plt.legend()
plt.show()

gates_output(gates_release)
