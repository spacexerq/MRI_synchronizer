# -*- coding: utf-8 -*-
"""
Created on Mon Jul  3 17:53:21 2023

@author: zilya
"""

import tkinter as tk
from math import ceil, floor 
import numpy as np
from types import SimpleNamespace
from datetime import datetime
import json
import pickle

def set_limits():
    
        # Задание общих аппаратных характкристик
    gamma = 42.576e6                   # Hz/T    Гиромагнитное отношение водорода
    G_amp_max_mT_m = 37.8                     # mT/m.   Максимальный градиент
    G_amp_max = G_amp_max_mT_m*1e-3*gamma   # Hz/m.   Максимальный градиент
    G_slew_max_T_m_s = 121                   # T/m/s.  Максимальная скорость нарастания
    G_slew_max = G_slew_max_T_m_s*gamma      # Hz/m/s. Максимальная скорость нарастания
    rf_raster_time = 1e-6              # s.      Растр РЧ импульса
    grad_raster_time = 10e-6           # s.      Растр градиентов
    tau_max = G_amp_max/G_slew_max     # s.      Максимальное время нарастания градиента с учетом макс скорости нарастания
    tau_max = np.ceil(tau_max/ grad_raster_time)*grad_raster_time
    tau = 250e-6                       # s.      Фиксированное время нарастания градиентов
    
    # Чтение заданных в интерфэйс значений 
    sl_nb = int(textBox1.get())
    sl_thkn = float(textBox2.get())
    sl_gap = float(textBox3.get())
    FoV_f = float(textBox4.get())
    FoV_ph_image = float(textBox5.get())
    Nf = int(textBox6.get())
    Np_image = int(textBox7.get())
    BW_pixel = float(textBox8.get())
    N_TE = int(textBox9.get())
    #TE = float(textBox9.get())
    TR = float(textBox10.get())
    alpha = float(textBox11.get())
    beta = float(textBox12.get())
    Conct = int(textBox13.get())
    ETL = int(textBox14.get())
    Average = int(textBox15.get())
    ph_over = float(textBox16.get())
    
    Np = Np_image*(1+ph_over/100)
    Np = round(Np)
    ph_over = (Np/Np_image-1)*100
    FoV_ph = FoV_ph_image*(1+ph_over/100)
    
    TR = np.ceil(TR/ grad_raster_time)*grad_raster_time
    
    
    # --- Параметры РЧ импульса -----
    rf_ringdown_time=20e-6
    rf_dead_time=100e-6
    adc_dead_time=10e-6
    
    # -  Библиотека time bandwidth products 
    t_BW_product_ex1 = 3.8
    t_BW_product_ref1 = 4.2
    t_BW_product_ex2 = 3.55
    t_BW_product_ref2 = 3.55
     
    # s. Библиотека длительностей импульсов
    # Длительности кратны растру градиента 10 мкс
    t_ex1 = 2.05e-3
    t_ref1 = 2.56e-3
    t_ex2 = 3.10e-3
    t_ref2 = 3.88e-3
    
    apodization = 0.27

    # ----- Подбор импульсов на основе выбранной толщины среза
    if sl_thkn > t_BW_product_ex1/(t_ex1*G_amp_max):  
       t_BW_product_ex = t_BW_product_ex1
       t_BW_product_ref = t_BW_product_ref1
    else:
       t_BW_product_ex = t_BW_product_ex2
       t_BW_product_ref = t_BW_product_ref2
        
    if sl_thkn > t_BW_product_ex2/(t_ex1*G_amp_max):
        t_ex = t_ex1
        t_ref = t_ref1 
    else:
        t_ex = t_ex2
        t_ref = t_ref2
     
    BW_ex_pulse = t_BW_product_ex/t_ex 
    
    # ---- Вычисление основных площадей и длительностей ------
    A_ex = BW_ex_pulse*(t_ex+tau)/sl_thkn
    A_read = Nf/FoV_f
    A_ph = Np/2/FoV_ph
    A_reph = 0.25*A_ex
    A_pre = 1.5*A_read
    A_crs = 0.75*A_ex
    A_crr = A_read
    A_sps = 4*A_ex
    
    t_reph_min = A_reph/G_amp_max + tau
    t_pre_min = A_pre/G_amp_max + tau
    t_pre_min = max(t_reph_min,t_pre_min)
    t_pre_min = np.ceil(t_pre_min / grad_raster_time)*grad_raster_time 
      
    t_crs_min = A_crs/G_amp_max + tau
    t_crr_min = A_crr/G_amp_max  + tau
    t_ph_min = A_ph/G_amp_max  + tau
    t_cr_min = max(t_crs_min,t_crr_min,t_ph_min)
    t_cr_min = np.ceil(t_cr_min / grad_raster_time)*grad_raster_time
    
    t_read = 1/BW_pixel
    t_read = np.ceil(t_read / grad_raster_time)*grad_raster_time 
    
    t_sps = A_sps/G_amp_max  + tau
    t_sps = np.ceil(t_sps / grad_raster_time)*grad_raster_time 
    
    
    ES_1 = t_ref + t_read + 2*t_cr_min + 4*tau 
    ES_2 = t_ref + t_ex + 2*t_cr_min + 2*t_pre_min + 4*tau 
    ES = max(ES_1,ES_2) 
    ES = np.ceil(ES/ grad_raster_time)*grad_raster_time
    
    t_cr = 0.5*(ES-t_ref-t_read) - 2*tau 
    t_cr = np.ceil(t_cr / grad_raster_time)*grad_raster_time
    
    t_pre = 0.5*(ES-t_ref-t_ex) - 2*tau 
    t_pre = np.ceil(t_pre / grad_raster_time)*grad_raster_time
    
    TE = N_TE*ES

    N_TE_min = 1
    N_TE_max = ETL

    TR_min = ES*(ETL) + t_ex/2 + t_read/2 + t_sps + 2*tau
    TR_max = 10
    TR_min = np.ceil(TR_min/ grad_raster_time)*grad_raster_time
    
    delay_TR = TR-TR_min
    
    Conc_min = ceil(sl_nb*TR_min/TR)
    Conc_max = sl_nb
    ETL_min = 1;
    ETL_max1 = 16;
    ETL_max2 = floor((TR-t_ex/2-t_read/2-t_sps-2*tau)/ES)
    ETL_max = min(ETL_max1,ETL_max2)
    

    Nf_min = 1
    #Nf_max1 = FoV_f*G_amp_max*(t_cr-tau)
    Nf_max = min(FoV_f*G_amp_max*(t_read),1600)
    #Nf_max = min(Nf_max1,Nf_max2)
    Np_min = 1
    Np_max1 = Nf
    Np_max2 = 2*FoV_ph*G_amp_max*(t_cr-tau)
    Np_max = min(Np_max1,Np_max2)

    FoV_min = 25e-3
    FoV_f_max = 450e-3
    FoV_f_min = Nf/(G_amp_max*t_read)
    #FoV_f_min2 = Nf/(G_amp_max*(t_cr-tau))
    FoV_f_min = max(FoV_f_min,FoV_min)
    FoV_p_min = FoV_f_min
    FoV_p_max = FoV_f
    
    Np_image_min = Np_min
    Np_image_max = min(Np_max, Nf_max/(1+ph_over/100))
    
    FoV_p_image_min = FoV_p_min/(1+ph_over/100)
    FoV_p_image_max = min(FoV_p_max,FoV_f_max/(1+ph_over/100))
    
    
   # BW_pixel_min = 1/(ES-t_ref-2*t_cr_min-4*tau) 
    # BW_pixel_min = max(ceil(BW_pixel_min), 40)
    BW_pixel_min = 40
    BW_pixel_max1 = 1780
    BW_pixel_max2 = FoV_f*G_amp_max/Nf
    BW_pixel_max = min (BW_pixel_max1,BW_pixel_max2)
    
    sl_thkn_min1 = t_BW_product_ex2/(t_ex2*G_amp_max)  
    sl_thkn_min2 = t_BW_product_ex2*(t_ex2+tau)/((t_cr-tau)*G_amp_max)  
    sl_thkn_min = max(sl_thkn_min1,sl_thkn_min2)
    sl_thkn_max = 20e-3
    sl_gap_min, sl_gap_max = 0, 800
    
    Average_min, Average_max = 1, 10
    
    ph_over_min, ph_over_max = 0, 100

        
    
    label2_1.configure(text = "min = " + str(ceil(sl_thkn_min*10000)/10000))
    label2_2.configure(text = "max = " + str(sl_thkn_max))
    label3_1.configure(text = "min = " + str(sl_gap_min))
    label3_2.configure(text = "max = " + str(sl_gap_max))
    label4_1.configure(text = "min = " + str(ceil(FoV_f_min*1000)/1000))
    label4_2.configure(text = "max = " + str(FoV_f_max))
    label5_1.configure(text = "min = " + str(ceil(FoV_p_image_min*1000)/1000))
    label5_2.configure(text = "max = " + str(ceil(FoV_p_image_max*1000)/1000))
    label6_1.configure(text = "min = " + str(Nf_min))
    label6_2.configure(text = "max = " + str(ceil(Nf_max)))
    label7_1.configure(text = "min = " + str(ceil(Np_image_min*1000)/1000))
    label7_2.configure(text = "max = " + str(ceil(Np_image_max*1000)/1000))
    label8_1.configure(text = "min = " + str(BW_pixel_min))
    label8_2.configure(text = "max = " + str(ceil(BW_pixel_max*1000)/1000))
    label9_1.configure(text = "min = " + str(N_TE_min))
    label9_2.configure(text = "max = " + str(N_TE_max))
    label10_1.configure(text = "min = " + str(ceil(TR_min*1000)/1000))
    label10_2.configure(text = "max = " + str(TR_max))
    label13_1.configure(text = "min = " + str(Conc_min))
    label13_2.configure(text = "max = " + str(Conc_max))
    label14_1.configure(text = "min = " + str(ETL_min))
    label14_2.configure(text = "max = " + str(ETL_max))
    label15_1.configure(text = "min = " + str(Average_min))
    label15_2.configure(text = "max = " + str(Average_max))
    label16_1.configure(text = "min = " + str(ph_over_min))
    label16_2.configure(text = "max = " + str(ph_over_max))
    
    global param
    param = SimpleNamespace()
    param.G_amp_max =G_amp_max_mT_m
    param.G_slew_max = G_slew_max_T_m_s
    param.gamma = gamma
    param.grad_raster_time = grad_raster_time
    param.rf_raster_time = rf_raster_time
    
    
    param.t_BW_product_ex = t_BW_product_ex
    param.t_BW_product_ref = t_BW_product_ref
    param.t_ex = t_ex
    param.t_ref = t_ref
    param.rf_ringdown_time = rf_ringdown_time
    param.rf_dead_time = rf_dead_time
    param.adc_dead_time = adc_dead_time
    param.aapodization = apodization
    
    param.dG = tau
    param.sl_nb = sl_nb
    param.sl_thkn = sl_thkn   
    param.sl_gap = sl_gap
    param.FoV_f = FoV_f
    param.FoV_ph = FoV_ph
    param.Nf = Nf
    param.Np = Np   
    param.BW_pixel = BW_pixel
    param.TE = TE
    param.N_TE = N_TE  #то на коком номере эхо истинное время эхо(в центр к-пр)
    param.ES = ES
    param.TR = TR
    param.FA = alpha
    param.RA = beta
    param.conct = Conct
    param.ETL = ETL
    param.Average = Average
    param.ph_over = ph_over
    
    param.delay_TR = delay_TR

   # tk.Label(win, text = 't = ' + str(round(time//60)) + ':'+ str(round(time%60)) + 'min').grid(row=13, column=2,sticky = "E")
    tk.Label(win, text = 'ES = ' + str(ceil(ES*10000)/10000) + ' s').grid(row=17, column=2,sticky = "E")
    tk.Label(win, text = 'TE = ' + str(ceil(TE*10000)/10000) + ' s').grid(row=18, column=2,sticky = "E")
    #tk.Label(win, text =  str(ETL_max2) + 's').grid(row=18, column=2,sticky = "E")

def save_param():
    
    output_filename = "test_1701" 
    # output_filename = "TSE_T1" + datetime.now().strftime("%Y%m%d_%H%M%S")
    file = open(output_filename + ".json", 'w')
    json.dump(param.__dict__, file, indent = 4)
    file.close() 
    
    # output_filename = "TSE_" + datetime.now().strftime("%Y%m%d_%H%M%S")
    # file = open(output_filename, 'wb')
    # pickle.dump(param, file)
    # file.close()

### Defoult values ###
# --------------------------------
sl_nb = 1          # number of slices                     
sl_thkn = 5e-3     # Slice thickness, m                   
sl_gap = 100       # Slice gap, %   
                      
#---------------------------------
FoV_f = 32e-3    # FoV in freq encoding direction, m                  
FoV_ph = 32e-3   # FoV in ph encoding direction, m               
Nf = 16          # Freq direction resolution                         
Np = 16          # Resolution in ph encoding direction              

#----------------------------------
BW_pixel = 500    # Pixel BW 
N_TE = 1          # number of echo in center of k-space
# TE = 20e-3         # Echo time, s
TR =  500e-3       # Repetition time, s
alpha = 90        # Flip angle, degree 
beta = 180
ETL = 8
Conct = 1
Average = 1
ph_over = 0

           

win = tk.Tk()
win.title('TSE')
win.geometry("350x410+100+100")
win.resizable(False,False)


# number of slices
sl_nb_min,sl_nb_max = 1, 103
tk.Label(win, text = 'Slices').grid(row=0, column=0,sticky = "E")
textBox1 = tk.Spinbox(from_=sl_nb_min, to=sl_nb_max, width = 4)
textBox1.grid(row=0, column=1)
label1_1 = tk.Label(win, text = "min = " + str(sl_nb_min))
label1_1.grid(row=0, column=2)
label1_2 = tk.Label(win, text = "max = " + str(sl_nb_max))
label1_2.grid(row=0, column=3)

# Slice thickness, m  
tk.Label(win, text = 'Slice thickness, m').grid(row=1, column=0,sticky = "E")
textBox2 = tk.Entry(width = 6)
textBox2.insert(0, sl_thkn)
textBox2.grid(row=1, column=1)
label2_1 = tk.Label(win, text = "min = " )
label2_1.grid(row=1, column=2)
label2_2 = tk.Label(win, text = "max = " )
label2_2.grid(row=1, column=3)

# Slice gap, %  
tk.Label(win, text = 'Slice gap, % ').grid(row=2, column=0,sticky = "E")
textBox3 = tk.Entry(width = 6)
textBox3.insert(0, sl_gap)
textBox3.grid(row=2, column=1)
label3_1 = tk.Label(win, text = "min = ")
label3_1.grid(row=2, column=2)
label3_2 = tk.Label(win, text = "max = ")
label3_2.grid(row=2, column=3)

# FoV in freq encoding direction, m
tk.Label(win, text = 'FoV read, m').grid(row=3, column=0,sticky = "E")
textBox4 = tk.Entry(width = 6)
textBox4.insert(0, FoV_f)
textBox4.grid(row=3, column=1)
label4_1 = tk.Label(win, text = "min = ")
label4_1.grid(row=3, column=2)
label4_2 = tk.Label(win, text = "max = ")
label4_2.grid(row=3, column=3)


# FoV in ph encoding direction, m  
label5 = tk.Label(win, text = 'FoV phase, m')
label5.grid(row=4, column=0,sticky = "E")
textBox5 = tk.Entry(width = 6)
textBox5.insert(0, FoV_f)
textBox5.grid(row=4, column=1)
label5_1 = tk.Label(win, text = "min = ")
label5_1.grid(row=4, column=2)
label5_2 = tk.Label(win, text = "max = ")
label5_2.grid(row=4, column=3)

# Freq direction resolution
label6 = tk.Label(win, text = 'Read resolution')
label6.grid(row=5, column=0,sticky = "E")
textBox6 = tk.Entry(width = 6)
textBox6.insert(0, Nf)
textBox6.grid(row=5, column=1)
label6_1 = tk.Label(win, text = "min = ")
label6_1.grid(row=5, column=2)
label6_2 = tk.Label(win, text = "max = ")
label6_2.grid(row=5, column=3)

# Ph direction resolution
label7 = tk.Label(win, text = 'Phase resolution')
label7.grid(row=6, column=0,sticky = "E")
textBox7 = tk.Entry(width = 6)
textBox7.insert(0, Np)
textBox7.grid(row=6, column=1)
label7_1 = tk.Label(win, text = "min = ")
label7_1.grid(row=6, column=2)
label7_2 = tk.Label(win, text = "max = ")
label7_2.grid(row=6, column=3)

# Pixel BW
label8 = tk.Label(win, text = 'Pixel BW, Hz')
label8.grid(row=7, column=0,sticky = "E")
textBox8 = tk.Entry(width = 6)
textBox8.insert(0, BW_pixel)
textBox8.grid(row=7, column=1)
label8_1 = tk.Label(win, text = "min = ")
label8_1.grid(row=7, column=2)
label8_2 = tk.Label(win, text = "max = ")
label8_2.grid(row=7, column=3)

# TE, s
label9 = tk.Label(win, text = 'N_TE')
label9.grid(row=8, column=0,sticky = "E")
textBox9 = tk.Entry(width = 6)
textBox9.insert(0, N_TE)
textBox9.grid(row=8, column=1)
label9_1 = tk.Label(win, text = "min = ")
label9_1.grid(row=8, column=2)
label9_2 = tk.Label(win, text = "max = ")
label9_2.grid(row=8, column=3)

# TR, s
label10 = tk.Label(win, text = 'TR, s')
label10.grid(row=9, column=0,sticky = "E")
textBox10 = tk.Entry(width = 6)
textBox10.insert(0, TR)
textBox10.grid(row=9, column=1)
label10_1 = tk.Label(win, text = "min = ")
label10_1.grid(row=9, column=2)
label10_2 = tk.Label(win, text = "max = ")
label10_2.grid(row=9, column=3)

# Flip angle, degree
label11= tk.Label(win, text = 'Flip angle, degree')
label11.grid(row=10, column=0,sticky = "E")
textBox11 = tk.Entry(width = 6)
textBox11.insert(0, alpha)
textBox11.grid(row=10, column=1)
angle_min, angle_max = 1, 90
label11_1 = tk.Label(win, text = "min = "+ str(angle_min))
label11_1.grid(row=10, column=2)
label11_2 = tk.Label(win, text = "max = "+ str(angle_max))
label11_2.grid(row=10, column=3)

# Refocusing angle, degree
label12= tk.Label(win, text = 'Refocusing angle, degree')
label12.grid(row=11, column=0,sticky = "E")
textBox12 = tk.Entry(width = 6)
textBox12.insert(0, beta)
textBox12.grid(row=11, column=1)
beta_min, beta_max = 110, 180
label12_1 = tk.Label(win, text = "min = "+ str(beta_min))
label12_1.grid(row=11, column=2)
label12_2 = tk.Label(win, text = "max = "+ str(beta_max))
label12_2.grid(row=11, column=3)

# Concatenations
label13= tk.Label(win, text = 'Concatenations')
label13.grid(row=12, column=0,sticky = "E")
textBox13 = tk.Entry(width = 6)
textBox13.insert(0, Conct)
textBox13.grid(row=12, column=1)
label13_1 = tk.Label(win, text = "min = ")
label13_1.grid(row=12, column=2)
label13_2 = tk.Label(win, text = "max = ")
label13_2.grid(row=12, column=3)

# ETL
label14= tk.Label(win, text = 'ETL')
label14.grid(row=13, column=0,sticky = "E")
textBox14 = tk.Entry(width = 6)
textBox14.insert(0, ETL)
textBox14.grid(row=13, column=1)
label14_1 = tk.Label(win, text = "min = ")
label14_1.grid(row=13, column=2)
label14_2 = tk.Label(win, text = "max = ")
label14_2.grid(row=13, column=3)

# Averages
label15= tk.Label(win, text = 'Averages')
label15.grid(row=14, column=0,sticky = "E")
textBox15 = tk.Entry(width = 6)
textBox15.insert(0, Average)
textBox15.grid(row=14, column=1)
label15_1 = tk.Label(win, text = "min = ")
label15_1.grid(row=14, column=2)
label15_2 = tk.Label(win, text = "max = ")
label15_2.grid(row=14, column=3)

# phase oversampling
label16= tk.Label(win, text = 'Ph oversampling,%')
label16.grid(row=15, column=0,sticky = "E")
textBox16 = tk.Entry(width = 6)
textBox16.insert(0, ph_over)
textBox16.grid(row=15, column=1)
label16_1 = tk.Label(win, text = "min = ")
label16_1.grid(row=15, column=2)
label16_2 = tk.Label(win, text = "max = ")
label16_2.grid(row=15, column=3)

set_limits()

btn1 = tk.Button(win, text = 'Set', command = set_limits)
btn1.grid(row=16, column=0)

btn1 = tk.Button(win, text = 'Save', command = save_param)
btn1.grid(row=16, column=3)



win.mainloop()


