#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 24 17:08:59 2022

@author: phillips
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

FONT_SIZE = 8

plt.rc('font', size=FONT_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=FONT_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=FONT_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=FONT_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=FONT_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=FONT_SIZE)    # legend fontsize
plt.rc('figure', titlesize=FONT_SIZE)  # fontsize of the figure title
plt.rcParams['pdf.fonttype'] = 42

import tensorflow_probability as tfp
import gpflow
import tensorflow as tf
from gpflow.utilities import print_summary
from gpflow.utilities import print_summary, set_trainable, to_default_float

from scipy.signal import find_peaks
from scipy.interpolate import CubicSpline

#%%

def load_data(file_name):
    df = pd.read_csv(file_name).fillna(0)  
    data_cols = [col for col in df if col.startswith('Cell')]
    bckgd_cols = [col for col in df if col.startswith('Background')]
    time = df['Time (h)'].values[:,None]

    bckgd = df[bckgd_cols].values
    M = np.shape(bckgd)[1]
    
    bckgd_length = np.zeros(M,dtype=np.int32)
    
    for i in range(M):
        bckgd_curr = bckgd[:,i]
        bckgd_length[i] = np.max(np.nonzero(bckgd_curr))
        
    y_all = df[data_cols].values  
    
    N = np.shape(y_all)[1]
    
    y_all = df[data_cols].values
    np.max(np.nonzero(y_all))
    
    y_length = np.zeros(N,dtype=np.int32)
    
    for i in range(N):
        y_curr = y_all[:,i]
        y_length[i] = np.max(np.nonzero(y_curr))    

    return time, bckgd, bckgd_length, M, y_all, y_length, N
    
#file_name = '/Users/phillips/Documents/GPosc python tutorial/Hes1_example.csv' #add your file name here   



# original Hes1 cells
# num = xlsread('Hes1ublucF5_sparse_combine_30minexp_forBayesian.xlsx',1);

def optimised_background_model(X,Y):
    
    k = gpflow.kernels.SquaredExponential()
    m = gpflow.models.GPR(data=(X, Y), kernel=k, mean_function=None) 
    m.kernel.lengthscales = gpflow.Parameter(to_default_float(7.1), transform=tfp.bijectors.Softplus(low=to_default_float(7.)))
    opt = gpflow.optimizers.Scipy()
    opt_logs = opt.minimize(m.training_loss, m.trainable_variables, options=dict(maxiter=100))

    return m

def background_std_vec(time, bckgd, bckgd_length, M):
    std_vec = np.zeros(M)

    fig = plt.figure(figsize=(15/2.54,15/2.54))
    
    for i in range(M):
        
        X = time[:bckgd_length[i]]
        Y = bckgd[:bckgd_length[i],i,None]  
        Y = Y - np.mean(Y)
        
        m = optimised_background_model(X,Y)
        
        mean, var = m.predict_y(X)
        
        plt.subplot(2,2,i+1)
        plt.plot(X,Y)
        plt.plot(X, mean,'k')
        plt.plot(X, mean + 2*var**0.5,'r')
        plt.plot(X, mean - 2*var**0.5,'r')
        
        if i%2==0:
            plt.ylabel('Luminescence (AU)')
        if i>=2:
            plt.xlabel('Time (h)')
        
        std_vec[i] = m.likelihood.variance**0.5
        
    return std_vec

def detrend_cell(X,Y,detrend_lengthscale):
    
    k_trend = gpflow.kernels.SquaredExponential()
    m = gpflow.models.GPR(data=(X, Y), kernel=k_trend, mean_function=None)

    m.kernel.lengthscales = gpflow.Parameter(to_default_float(detrend_lengthscale+0.1), transform=tfp.bijectors.Softplus(low=to_default_float(detrend_lengthscale)))
    
    opt = gpflow.optimizers.Scipy()
    opt_logs = opt.minimize(m.training_loss, m.trainable_variables, options=dict(maxiter=100))

    mean, var = m.predict_f(X)
    
    Y_detrended = Y - mean
    Y_detrended = Y_detrended-np.mean(Y_detrended)
    
    return k_trend, mean, var, Y_detrended

def fit_models(X,Y,noise,K):
    
    OU_LL_list, OU_param_list, OUosc_LL_list, OUosc_param_list = [[] for _ in range(4)]
    
    for k in range(K):
        
        try:
        
            k_ou = gpflow.kernels.Matern12()
        
            m = gpflow.models.GPR(data=(X, Y), kernel=k_ou, mean_function=None)
            m.kernel.variance.assign(np.random.uniform(0.1,2.0))
            m.kernel.lengthscales.assign(np.random.uniform(0.1,2.0))
            m.likelihood.variance.assign(noise**2)
            gpflow.set_trainable(m.likelihood.variance,False)  
            opt = gpflow.optimizers.Scipy()
            opt_logs = opt.minimize(m.training_loss, m.trainable_variables, options=dict(maxiter=100))
            
            nlmlOU = m.log_posterior_density()
            
            OU_LL = nlmlOU
            OU_LL_list.append(OU_LL)
            OU_param_list.append(k_ou)
            
            k_ou_osc = gpflow.kernels.Matern12()*gpflow.kernels.Cosine()
        
            m = gpflow.models.GPR(data=(X, Y), kernel=k_ou_osc, mean_function=None)
            m.likelihood.variance.assign(noise**2)
            gpflow.set_trainable(m.likelihood.variance,False)  
            gpflow.set_trainable(m.kernel.kernels[1].variance,False)  
            m.kernel.kernels[0].variance.assign(np.random.uniform(0.1,2.0))
            m.kernel.kernels[0].lengthscales.assign(np.random.uniform(0.1,2.0))
            m.kernel.kernels[1].lengthscales.assign(np.random.uniform(0.1,4.0))
            opt = gpflow.optimizers.Scipy()
            opt_logs = opt.minimize(m.training_loss, m.trainable_variables, options=dict(maxiter=100))
            
            nlmlOSC = m.log_posterior_density()#opt_logs.fun
            
            OU_osc_LL = nlmlOSC
            OUosc_LL_list.append(OU_osc_LL)
            OUosc_param_list.append(k_ou_osc)
            
        except:
            pass
         
    LLR = 100*2*(np.max(OUosc_LL_list)-np.max(OU_LL_list))/len(Y)
    BIC_OUosc = -2*np.max(OUosc_LL_list)+3*np.log(len(Y))
    BIC_OU = -2*np.max(OU_LL_list)+2*np.log(len(Y))
    BICdiff = BIC_OU-BIC_OUosc
    k_ou = OU_param_list[np.argmax(OU_LL_list)]
    k_ou_osc = OUosc_param_list[np.argmax(OUosc_LL_list)]
    
    cov_ou_osc = OUosc_param_list[0](X).numpy()[0,:]
    peaks, _ = find_peaks(cov_ou_osc, height=0)
    
    if len(peaks)!=0:
        period = X[peaks[0]]
    else:
        period = 0
    
    return LLR, BICdiff, k_ou, k_ou_osc, period

def plot_model_fits(cell, x_curr,y_curr,mean_trend,noise,LLR, k_trend, k_ou, k_ou_osc, period):
    fig = plt.figure(figsize=(12/2.54,8/2.54))
    plt.plot(x_curr,y_curr)
    plt.plot(x_curr,mean_trend,'k--',alpha=0.5)
    plt.xlabel('Time (hours)')
    plt.ylabel('Luminescence (normalised) (AU)')
    plt.title('Cell '+str(cell)+' , LLR = '+f'{LLR:.1f}')
    

def analyse_cells(time,y_all, y_length, N):
    
    noise_list, detrend_param_list, LLR_list, BICdiff_list, OU_param_list, OUosc_param_list, period_list = [[] for _ in range(7)]
    
    for cell in range(N):
    
        print(cell)
        
        x_curr = time[:y_length[cell]]    
        y_curr = y_all[:y_length[cell],cell,None]
        noise = std/np.std(y_curr)
        y_curr=(y_curr-np.mean(y_curr))/np.std(y_curr)
    
        k_trend, mean_trend, var_trend, Y_detrended = detrend_cell(x_curr,y_curr,7.0)
        
        LLR, BICdiff, k_ou, k_ou_osc, period = fit_models(x_curr,Y_detrended,noise,10)
        
        if cell==0:
            plot_model_fits(cell, x_curr,y_curr,mean_trend,noise,LLR, k_trend, k_ou, k_ou_osc, period) 
       
        noise_list.append(noise)
        detrend_param_list.append(k_trend)
        LLR_list.append(LLR)
        BICdiff_list.append(BICdiff)
        OU_param_list.append(k_ou)
        OUosc_param_list.append(k_ou_osc)
        period_list.append(period) 
    
    return noise_list, detrend_param_list, LLR_list, BICdiff_list, OU_param_list, OUosc_param_list, period_list

file_name = '/Users/phillips/Documents/GPosc python tutorial/Hes1_batch1.csv' #add your file name here  
time, bckgd, bckgd_length, M, y_all, y_length, N = load_data(file_name)

std_vec = background_std_vec(time, bckgd, bckgd_length, M)
std = np.mean(std_vec)  # the estimated standard deviation of the experimental noise, averaged over all background

noise_list_HES_batch1, detrend_param_list_HES_batch1, LLR_list_HES_batch1, BICdiff_list_HES_batch1, OU_param_list_HES_batch1, OUosc_param_list_HES_batch1, period_list_HES_batch1 = analyse_cells(time,y_all, y_length, N)
N_batch1 = N
y_length_batch1 = y_length

file_name = '/Users/phillips/Documents/GPosc python tutorial/Hes1_batch2.csv' #add your file name here  
time, bckgd, bckgd_length, M, y_all, y_length, N = load_data(file_name)

std_vec = background_std_vec(time, bckgd, bckgd_length, M)
std = np.mean(std_vec)  # the estimated standard deviation of the experimental noise, averaged over all background

noise_list_HES_batch2, detrend_param_list_HES_batch2, LLR_list_HES_batch2, BICdiff_list_HES_batch2, OU_param_list_HES_batch2, OUosc_param_list_HES_batch2, period_list_HES_batch2 = analyse_cells(time,y_all, y_length, N)
N_batch2 = N
y_length_batch2 = y_length


N_HES_tot = N_batch1+N_batch2
y_length_HES_tot = np.concatenate((y_length_batch1,y_length_batch2),0)
noise_list_HES_tot = noise_list_HES_batch1+noise_list_HES_batch2
detrend_param_list_HES_tot = detrend_param_list_HES_batch1+detrend_param_list_HES_batch2
OU_param_list_HES_tot = OU_param_list_HES_batch1+OU_param_list_HES_batch2
LLR_list_HES_tot = LLR_list_HES_batch1+LLR_list_HES_batch2



file_name = '/Users/phillips/Documents/GPosc python tutorial/Control_batch1.csv' #add your file name here  
time, bckgd, bckgd_length, M, y_all, y_length, N = load_data(file_name)

std_vec = background_std_vec(time, bckgd, bckgd_length, M)
std = np.mean(std_vec)  # the estimated standard deviation of the experimental noise, averaged over all background

noise_list_Control_batch1, detrend_param_list_Control_batch1, LLR_list_Control_batch1, BICdiff_list_Control_batch1, OU_param_list_Control_batch1, OUosc_param_list_Control_batch1, period_list_Control_batch1 = analyse_cells(time,y_all, y_length, N)
N_batch1 = N
y_length_batch1 = y_length

file_name = '/Users/phillips/Documents/GPosc python tutorial/Control_batch2.csv' #add your file name here  
time, bckgd, bckgd_length, M, y_all, y_length, N = load_data(file_name)

std_vec = background_std_vec(time, bckgd, bckgd_length, M)
std = np.mean(std_vec)  # the estimated standard deviation of the experimental noise, averaged over all background

noise_list_Control_batch2, detrend_param_list_Control_batch2, LLR_list_Control_batch2, BICdiff_list_Control_batch2, OU_param_list_Control_batch2, OUosc_param_list_Control_batch2, period_list_Control_batch2 = analyse_cells(time,y_all, y_length, N)
N_batch2 = N
y_length_batch2 = y_length



N_Control_tot = N_batch1+N_batch2
y_length_Control_tot = np.concatenate((y_length_batch1,y_length_batch2),0)
noise_list_Control_tot = noise_list_Control_batch1+noise_list_Control_batch2
detrend_param_list_Control_tot = detrend_param_list_Control_batch1+detrend_param_list_Control_batch2
OU_param_list_Control_tot = OU_param_list_Control_batch1+OU_param_list_Control_batch2
LLR_list_Control_tot = LLR_list_Control_batch1+LLR_list_Control_batch2



N_tot = N_HES_tot+N_Control_tot
y_length_tot = np.concatenate((y_length_HES_tot,y_length_Control_tot),0)
noise_list_tot = noise_list_HES_tot+noise_list_Control_tot
detrend_param_list_tot = detrend_param_list_HES_tot+detrend_param_list_Control_tot
OU_param_list_tot = OU_param_list_HES_tot+OU_param_list_Control_tot

#%%

def synthetic_cell_LLRs(N,y_length,noise_list,detrend_param_list,OU_param_list):
    repeats = 10

    LLR_list_synth = []
    
    for cell in range(N):
        
        print("Progress: {0}/{1}".format(cell,N))
        
        X = time[:y_length[cell]]  
        noise = noise_list[cell]
        
        k_se = detrend_param_list[cell]
        k_ou = OU_param_list[cell]
        k_white = gpflow.kernels.White(variance = noise**2)
        
        k_synth = k_se + k_ou + k_white
        
        for repeat in range(repeats):
        
            y_synth = np.random.multivariate_normal(np.zeros(len(X)), k_synth(X)).reshape(-1,1)        
            k_trend, mean_trend, var_trend, Y_detrended = detrend_cell(X,y_synth,7.0)
            LLR, BICdiff, k_ou, k_ou_osc, period = fit_models(X,Y_detrended,noise,10)
            LLR_list_synth.append(LLR)
    return LLR_list_synth        

LLR_list_synth = synthetic_cell_LLRs(N_tot,y_length_tot,noise_list_tot,detrend_param_list_tot,OU_param_list_tot)

#%%

def calc_q_vals(LLR_list,LLR_list_synth):
    LLR_array = np.array(LLR_list)
    LLR_synth_array = np.array(LLR_list_synth)
    
    # LLRs can be tiny and just negative - this just sets them to zero
    LLR_array[LLR_array<0] = 0
    LLR_synth_array[LLR_synth_array<0] = 0
    
    LLR_combined = np.concatenate((LLR_array,LLR_synth_array),0)
    
    upper = np.max(LLR_combined)
    lower1 = np.min(LLR_combined)
    lower = upper - 0.9*(upper-lower1);
    
    grid = np.linspace(lower,upper,20)
    
    piest = np.zeros_like(grid)
    
    for i, cutoff in enumerate(grid):
        num = sum(LLR_array<cutoff)/len(LLR_array)
        denom = sum(LLR_synth_array<cutoff)/len(LLR_synth_array)
        piest[i] =  num/denom
        
    xx = np.linspace(lower,upper,100);
       
    cs = CubicSpline(grid, piest)
    yy = cs(xx)
    
    piGUESS1 = yy[0]
    
    I = np.argsort(LLR_array)
    
    LLR_array_sorted = LLR_array[I]
    
    q1 = np.zeros_like(LLR_array_sorted)
    
    for i, thresh in enumerate(LLR_array_sorted):
        q1[i] = piGUESS1*(sum(LLR_synth_array>=thresh)/len(LLR_synth_array))/(sum(LLR_array_sorted>=thresh)/len(LLR_array_sorted))
        
    q_vals = q1[np.argsort(I)]  
    osc_filt = q_vals<0.05
    
    return q_vals, osc_filt

q_vals_HES1, osc_filt_HES1 = calc_q_vals(LLR_list_HES_tot,LLR_list_synth)
print("Number of Hes1 cells counted as oscillatory: {0}/{1}".format(sum(osc_filt_HES1),len(osc_filt_HES1)))

q_vals_Control, osc_filt_Control = calc_q_vals(LLR_list_Control_tot,LLR_list_synth)
print("Number of control cells counted as oscillatory: {0}/{1}".format(sum(osc_filt_Control),len(osc_filt_Control)))

#%%

fig = plt.figure(figsize=(18/2.54,6/2.54))

plt.subplot(1,3,1)
plt.hist(LLR_list_HES_tot,bins=np.linspace(0,40,40))
plt.xlabel('LLR')
plt.ylabel('Frequency')
plt.title('LLRs of HES1 cells')

plt.subplot(1,3,2)
plt.hist(LLR_list_Control_tot,bins=np.linspace(0,40,40))
plt.xlabel('LLR')
plt.ylabel('Frequency')
plt.title('LLRs of control cells')

plt.subplot(1,3,3)
plt.hist(LLR_list_synth,bins=np.linspace(0,40,40))
plt.xlabel('LLR')
plt.ylabel('Frequency')
plt.title('LLRs of synthetic non-oscillatory OU cells')

plt.tight_layout()
