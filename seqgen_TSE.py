#---------------------------------------------------------------------
# imports of the libraries
#---------------------------------------------------------------------
from math import pi
import numpy as np
import math
import json as j

from pypulseq.Sequence.sequence import Sequence
from pypulseq.calc_rf_center import calc_rf_center
from pypulseq.calc_duration import calc_duration
from pypulseq.make_adc import make_adc
from pypulseq.make_delay import make_delay
from pypulseq.make_sinc_pulse import make_sinc_pulse
from pypulseq.make_trapezoid import make_trapezoid
from pypulseq.make_extended_trapezoid import make_extended_trapezoid
from pypulseq.opts import Opts
from pypulseq.align import align
from pypulseq.traj_to_grad import traj_to_grad

from utilities import phase_grad_utils as pgu

#from seq2xml_fixed_delay import seq2xml


def TSE_k_space_fill(n_ex, ETL, k_steps, TE_eff_number):
    # function defines phase encoding steps for k space filling in liner order
    # with shifting according to the TE effective number

    global k_space_order_filing, k_space_list_with_zero, k_space_list_temp
    k_space_list_with_zero = []
    for i in range(ETL):
        k_space_list_with_zero.append((ETL - 1) * n_ex - i * n_ex)
    # print(k_space_list_with_zero)
    central_num = np.int32(k_steps / 2)
    # print(central_num)
    index_central_line = k_space_list_with_zero.index(central_num)
    shift = index_central_line - TE_eff_number + 1

    if shift > 0:
        a = k_space_list_with_zero[:shift]
        b = k_space_list_with_zero[shift:]
        k_space_list_with_zero = b + a
    elif shift < 0:
        a = k_space_list_with_zero[:shift]
        b = k_space_list_with_zero[shift:]
        k_space_list_with_zero = b + a

    k_space_order_filing = [k_space_list_with_zero]
    for i in range(n_ex - 1):
        k_space_list_temp = []
        for k in k_space_list_with_zero:
            k_space_list_temp.append(k + i + 1)
        k_space_order_filing.append(k_space_list_temp)

    return k_space_order_filing


def main(plot: bool, write_seq: bool, weightning):

    # Reading json file according to the weightning of the image
    if weightning == 'T1': #TODO: create general path
        with open('test_1701.json', 'rb') as f:
            params = j.load(f)

    elif weightning == 'T2':
        with open('C:\MRI_seq_files_mess\TSE_T2.json', 'rb') as f:
            params = j.load(f)

    elif weightning == 'T1':
        with open('C:\MRI_seq_files_mess\TSE_T1.json', 'rb') as f:
            params = j.load(f)
    else:
        print('Please choose image weightning')

    readout_time = round(1 / params['BW_pixel'], 8)

    # --------------------------
    # Set system limits
    # --------------------------

    scanner_parameters = Opts(
        max_grad=37.8,
        grad_unit="mT/m",
        max_slew=121,
        slew_unit="T/m/s",
        rf_ringdown_time=params['rf_ringdown_time'],
        rf_dead_time=params['rf_dead_time'],
        adc_dead_time=params['adc_dead_time'],
        rf_raster_time=params['rf_raster_time'],
        grad_raster_time=params['grad_raster_time'],
        block_duration_raster=params['grad_raster_time'],
        adc_raster_time=1 / (params['BW_pixel'] * params['Nf'])
    )
    seq = Sequence(scanner_parameters)

    #--------------------------
    # RF & Gradients
    #--------------------------

    rf90_phase = np.pi / 2
    rf180_phase = 0

    flip90 = round(params['FA'] * pi / 180, 3)
    flip180 = round(180 * pi / 180)
    rf90, gz_ex, _ = make_sinc_pulse(flip_angle=flip90, system=scanner_parameters, duration=params['t_ex'],
                                      slice_thickness=params['sl_thkn'], apodization=0.3,
                                      time_bw_product=round(params['t_BW_product_ex'], 8), return_gz=True)

    rf180, gz_ref, _ = make_sinc_pulse(flip_angle=flip180, system=scanner_parameters, duration=params['t_ref'],
                                        slice_thickness=params['sl_thkn'], apodization=0.3,
                                        time_bw_product=round(params['t_BW_product_ref'], 8), phase_offset=90 * pi / 180,
                                        return_gz=True)
    


    # Prepare RF offsets. This is required for multi-slice acquisition
    pulse_offsets90 = (np.linspace(0.0, params['sl_nb'] - 1.0, np.int16(params['sl_nb'])) - 0.5 * (
                np.double(params['sl_nb']) - 1.0)) * (
                                  params['sl_thkn'] * (100.0 + params['sl_gap']) / 100.0) * gz_ex.amplitude
    pulse_offsets180 = (np.linspace(0.0, params['sl_nb'] - 1.0, np.int16(params['sl_nb'])) - 0.5 * (
                np.double(params['sl_nb']) - 1.0)) * (
                                   params['sl_thkn'] * (100.0 + params['sl_gap']) / 100.0) * gz_ref.amplitude

    # slice selective gradient drafts for complex gradient blocks
    t_exwd = params['t_ex'] + scanner_parameters.rf_ringdown_time + scanner_parameters.rf_dead_time
    t_refwd = params['t_ref'] + scanner_parameters.rf_ringdown_time + scanner_parameters.rf_dead_time

    gz90 = make_trapezoid(channel="z", system=scanner_parameters, amplitude=gz_ex.amplitude,
                          flat_time=t_exwd, rise_time=params['dG'])
    gz180 = make_trapezoid(channel="z", system=scanner_parameters, amplitude=gz_ref.amplitude,
                           flat_time=t_refwd, rise_time=params['dG'])

    # generate basic gx readout gradient - G_read
    k_read = np.double(params['Nf']) / np.double(params['FoV_f'])
    t_gx = np.ceil(readout_time / scanner_parameters.grad_raster_time) * scanner_parameters.grad_raster_time
    gx = make_trapezoid(channel='x', system=scanner_parameters, flat_area=k_read,
                        flat_time=t_gx + 2 * scanner_parameters.adc_dead_time)

    # generate gx spoiler gradient - G_crr
    gx_spoil = make_trapezoid(channel='x', system=scanner_parameters, area=gx.area, flat_time=params['dG'],
                              rise_time=params['dG'])

    # read prephase gradient - G_pre
    gx_pre = make_trapezoid(channel="x", system=scanner_parameters, area=gx.area * 1.50,
                            rise_time=params['dG'])

    # rephase gradient draft after 90 RF pulse  - G_reph
    gz_reph = make_trapezoid(channel="z", system=scanner_parameters, area=gz_ex.area * 0.25,
                             flat_time=calc_duration(gx_pre), rise_time=params['dG'])

    # spoil gradient around 180 RF pulse - G_crs
    t_gz_spoil = np.ceil(
        params['t_ref'] / 2 / scanner_parameters.grad_raster_time) * scanner_parameters.grad_raster_time
    gz_spoil = make_trapezoid(channel='z', system=scanner_parameters, area=gz90.area * 0.75, rise_time=params['dG'],
                              flat_time=params['dG'])

    # spoil gradient G_sps
    gz_cr = make_trapezoid(channel='z', system=scanner_parameters, area=gz90.area * 4, rise_time=params['dG'])

    # Creation of ADC
    adc = make_adc(num_samples=params['Nf'], duration=t_gx, delay=scanner_parameters.adc_dead_time,
                   system=scanner_parameters)

    #--------------------------
    # k-space filling quantification
    #--------------------------

    k_phase = np.double(params['Np']) / np.double(params['FoV_ph'])
    k_steps_PE = pgu.create_k_steps(k_phase, np.int16(params['Np']))  # list of phase encoding gradients

    n_ex = math.floor(params['Np'] / params['ETL'])  # number of excitations
    #k_space_order_filing = TSE_k_space_fill(n_ex, np.int32(params['ETL']), np.int32(params['Np']), np.int32(
    #    params['N_TE']))  # TODO: to create additiolal functions on different k space order filling
    k_space_order_filing = TSE_k_space_fill(4, np.int32(8), np.int32(32), np.int32(1))
    k_space_order_filing

    #--------------------------
    # DELAYS
    #--------------------------

    block_duration = 0
    block_duration = max(calc_duration(rf90), calc_duration(gz90)) / 2 
    block_duration += max(calc_duration(gx_pre), calc_duration(gz_spoil))
    for i in range(np.int32(params['ETL']) - 1):
        block_duration += max(calc_duration(rf180), calc_duration(gz180))
        block_duration += max(calc_duration(gx_spoil), calc_duration(gz_spoil))
        block_duration += max(calc_duration(gx_spoil), calc_duration(adc))
        block_duration += calc_duration(gz_spoil)
    block_duration += max(calc_duration(rf180), calc_duration(gz_spoil))
    block_duration += max(calc_duration(gx_spoil), calc_duration(gz_spoil))
    block_duration += max(calc_duration(gx_spoil), calc_duration(adc))
    block_duration += calc_duration(gz_cr)

    #--------------------------
    # CONSTRUCT CONCATINATIONS timings
    #--------------------------

    # Quantification of Effective TE loop
    # eff_time = TE + calc_duration(gx) / 2 + max(calc_duration(gy_pre),calc_duration(gz_spoil)) + calc_duration(gx_spoil) + calc_duration(gz90) / 2
    eff_time = block_duration  # equal to previous!

    max_slices_per_TR = np.floor(params['TR'] / eff_time)
    required_concats = np.int32(np.ceil(params['sl_nb'] / max_slices_per_TR))
    slice_list = list(range(np.int32(params['sl_nb'])))
    slice_list = [slice_list[x::required_concats] for x in range(required_concats)]

    # Calculate the TR fillers
    tr_pauses = [(params['TR'] / np.double(len(x))) - eff_time for x in slice_list]
    tr_pauses = [
        max(seq.grad_raster_time, seq.rf_raster_time) * np.floor(x / max(seq.grad_raster_time, seq.rf_raster_time)) for
        x in tr_pauses]

    # Generate the TR fillers
    tr_fillers = [make_delay(x) for x in tr_pauses]

    # --------------------------
    # CONSTRUCT SEQUENCE
    # --------------------------

    for k in range(int(params['Average'])):  # Averages
        for curr_concat in range(required_concats):
            for phase_steps in k_space_order_filing:  # in stead of phase steps list of phase steps
                for curr_slice in range(np.int32(params['sl_nb'])):  # Slices
                    # Apply RF offsets
                    n_echo_temp = 0
                    rf90.freq_offset = pulse_offsets90[curr_slice]
                    rf180.freq_offset = pulse_offsets180[curr_slice]
                    # rf90.phase_offset = (rf90_phase - 2 * np.pi * rf90.freq_offset * calc_rf_center(rf90)[0])
                    # rf180.phase_offset = (rf180_phase - 2 * np.pi * rf180.freq_offset * calc_rf_center(rf180)[0])
                    print('curr_concat_' + str(curr_concat))
                    print('curr_slice_' + str(curr_slice))

                    seq.add_block(gz90, rf90)
                    seq.add_block(gz_reph, gx_pre)
                    for phase_step in phase_steps:
                        print('phase step_' + str(phase_step))
                        seq.add_block(gz180, rf180)
                        gy_pre = make_trapezoid(channel='y', system=scanner_parameters,
                                                area=-k_steps_PE[phase_step], duration=calc_duration(gz_spoil),
                                                rise_time=params['dG'])
                        print(k_steps_PE[phase_step])

                        seq.add_block(gz_spoil, gx_spoil, gy_pre)
                        seq.add_block(gx, adc)
                        n_echo_temp += 1
                        gy_pre = make_trapezoid(channel='y', system=scanner_parameters,
                                                area=k_steps_PE[phase_step], duration=calc_duration(gz_spoil),
                                                rise_time=params['dG'])
                        if n_echo_temp == np.int32(params['ETL']):
                            seq.add_block(gz_cr, gx_spoil, gy_pre)
                        else:
                            seq.add_block(gz_spoil, gx_spoil, gy_pre)
                    seq.add_block(tr_fillers[curr_concat])

    ok, error_report = seq.check_timing()
    if ok:
        print("Timing check passed successfully")
    else:
        print("Timing check failed. Error listing follows:")
        [print(e) for e in error_report]

    # ======
    # VISUALIZATION
    # ======
    if plot:
        seq.plot(time_range = [0, 0.02])
    # =========
    # WRITE .SEQ
    # =========
    if write_seq: #TODO: create general path
        if weightning == 'T1':
            seq.write('test_1701.seq')  # Save to disk
            #seq2xml(seq, seq_name='t1_TSE_matrx16x16_myGrad', out_folder='C:\\MRI_seq\\new_MRI_pulse_seq\\t1_TSE')

        elif weightning == 'T2':
            seq.write('TSE 16-16.seq')  # Save to disk
           # seq2xml(seq, seq_name='t2_TSE_matrx16x16_myGrad', out_folder='C:\\MRI_seq\\new_MRI_pulse_seq\\t2_TSE')

        elif weightning == 'PD':
            seq.write('TSE 16-16.seq')  # Save to disk
           # seq2xml(seq, seq_name='t2_TSE_matrx16x16_myGrad', out_folder='C:\\MRI_seq\\new_MRI_pulse_seq\\t2_TSE')

        else:
            print('Please choose image weightning')


if __name__ == "__main__":
    main(plot=True, write_seq=True, weightning='T1')
