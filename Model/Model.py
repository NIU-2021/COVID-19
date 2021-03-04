#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')
get_ipython().system('pip install mpld3')
import mpld3
mpld3.enable_notebook()


# In[2]:


train_data = pd.read_csv(import the file location address here (Files are available for each state (csv format)))
train_data 


# In[3]:


Last_Date = train_data .iloc[-1]["Date"]



# In[4]:


get_ipython().system('pip install lmfit')
from lmfit import Model


# In[5]:


def New_R_0(t, R_0_start, k, x0, R_0_end):
    return (R_0_start-R_0_end) / (1 + np.exp(-k*(-t+x0))) + R_0_end


# In[6]:


from scipy.integrate import odeint 
def Eqs_SEIR(y, t, N, beta, gamma, sigma):
    S, E, I, R = y
    dSdt = -beta(t) * S * I / N  
    dEdt = beta(t) * S * I / N - sigma * E 
    dIdt = sigma * E - gamma * I  
    dRdt = gamma * I  
    return dSdt, dEdt, dIdt, dRdt


# In[7]:


def Modified_SEIR(N, D, Max_d, CFR_O, CFR_S, R_0, **R_0_kwargs):
    '''
    R_0: callable
    '''
    # Initial value
    I0, R0, E0 = 0, 0, 1  
    S0 = N - I0 - R0 - E0

    gamma = 1.0 / D
    sigma = 1.0 / 5.0 

    # Defining contact rate
    def beta(t):
        return R_0(t, **R_0_kwargs) * gamma
    t = np.linspace(0, Max_d, Max_d)

    # Initial vector
    y0 = S0, E0, I0, R0
    ret = odeint(Eqs_SEIR, y0, t, args=(N, beta, gamma, sigma))
    S, E, I, R = ret.T

    def CFR(t):
        if t < 7:
            return CFR_O
        else:
            return CFR_O + CFR_S * (I[t - 7] / N)  
                
    # Inorporating deaths
    X = np.zeros(Max_d)
    for day in range(16, Max_d): 
        for valid_day in range(day-16):
            X[day] += CFR(valid_day) * beta(valid_day) * I[valid_day]
             
    return t, S, E, I, R, X, [R_0(t, **R_0_kwargs) for t in range(Max_d)], N, [CFR(t) for t in range(Max_d)]


# In[8]:


def plot_Modified_SEIR(t, S, E, I, R, X, R_0, N, CFR, x_ticks=None):
    # Plot the data   
    f, ax = plt.subplots(1,1,figsize=(10,4))
    #ax.plot(t, S, 'b', alpha=0.7, linewidth=2.5, label='Susceptible')
    ax.plot(t, E, 'y--', alpha=0.7, linewidth=2.5, label='Exposed')
    ax.plot(t, I, 'y', alpha=0.7, linewidth=2.5, label='Infected')
    ax.plot(t, X, 'r', alpha=0.7, linewidth=2.5, label='Dead')
    ax.plot(t, R, 'g', alpha=0.7, linewidth=2.5, label='Recovered')

    ax.set_xlabel('Time (days)')
    ax.title.set_text('SEIR-Model')
    ax.yaxis.set_tick_params(length=0)
    ax.xaxis.set_tick_params(length=0)

    if x_ticks is not None:
        ax.set_xticks(t[::21])
        ax.set_xticklabels(x_ticks[::21])    

    legend = ax.legend()
    legend.get_frame().set_alpha(0.8)
    for spine in ('top', 'right', 'bottom', 'left'):
        ax.spines[spine].set_visible(False)
    plt.savefig('SEIR-Model_NJ.png', dpi=600)    
    plt.show();
    
    f = plt.figure(figsize=(10,4))

    ax1 = f.add_subplot(121)
    ax1.plot(t, R_0, 'b--', alpha=0.7, linewidth=2.5, label='R_0')
 
    ax1.set_xlabel('Time (days)')
    ax1.title.set_text('R_0 over time')
    ax1.yaxis.set_tick_params(length=0)
    ax1.xaxis.set_tick_params(length=0)
    if x_ticks is not None:
        ax1.set_xticks(t[::35])
        ax1.set_xticklabels(x_ticks[::35])    
    legend = ax1.legend()
    legend.get_frame().set_alpha(0.8)
    for spine in ('top', 'right', 'bottom', 'left'):
        ax.spines[spine].set_visible(False)

    ax2 = f.add_subplot(122)
    ax2.plot(t, CFR, 'r--', alpha=0.7, linewidth=2.5, label='CFR')
    
    ax2.set_xlabel('Time (days)')
    ax2.title.set_text('CFR over time')
    ax2.yaxis.set_tick_params(length=0)
    ax2.xaxis.set_tick_params(length=0)
    if x_ticks is not None:
        ax2.set_xticks(t[::70])
        ax2.set_xticklabels(x_ticks[::70])
    legend = ax2.legend()
    legend.get_frame().set_alpha(0.5)
    for spine in ('top', 'right', 'bottom', 'left'):
        ax.spines[spine].set_visible(False)
    plt.savefig('R_0 and CFR_NJ.png', dpi=600)    

    plt.show();


# In[9]:


def plot_Modified_SEIR_1(t, S, E, I, R, X, R_0, N, CFR, x_ticks=None):
    # Plot the data on three separate curves for S(t), I(t) and R(t)
    f, ax = plt.subplots(1,1,figsize=(10,4))
    ax.plot(t, I, 'r', alpha=0.7, linewidth=4, label='Infected')
 
    ax.set_xlabel('Time (days)')
    ax.title.set_text('SEIR-Model with varying R_0 and CFR')

    ax.yaxis.set_tick_params(length=0)
    ax.xaxis.set_tick_params(length=0)

    if x_ticks is not None:
        ax.set_xticks(t[::21])
        ax.set_xticklabels(x_ticks[::21])    

    legend = ax.legend()
    legend.get_frame().set_alpha(0.5)
    for spine in ('top', 'right', 'bottom', 'left'):
        ax.spines[spine].set_visible(False)
    plt.savefig('Infected vs time.png', dpi=600)
    plt.show();


# In[10]:


def Fit_Modified_SEIR(CN, missing_days=0, RN=None, fit_method="least_squares", **R_0_kwargs):

    if RN is not None:
        y_data = train_data [(train_data ["State"] == CN)].Fatalities.values
    else:
        if len(train_data ["State"] == CN) > len(train_data ["State"] == "United States"):  
            y_data = train_data [(train_data ["State"] == CN)].Fatalities.values
        else:
            y_data = train_data [train_data ["State"] == CN].Fatalities.values
        
    Max_d = len(train_data.groupby("Date").sum().index) + missing_days 
    y_data = np.concatenate((np.zeros(missing_days), y_data))
    
    # State Population
    N = "please import the state population here"
    
    # Infectious period (please read the manuscript for more information)
    D = 10

    x_data = np.linspace(0, Max_d - 1, Max_d, dtype=int)

    # Solver
    def Modified_SEIR_deaths(x, N, D, CFR_O, CFR_S, R_0_delta, **R_0_kwargs):
        t_, S_, E_, I_, R_, X, R_0_, N_, CFR_ = Modified_SEIR(N, D, Max_d, CFR_O, CFR_S, R_0=New_R_0, **R_0_kwargs)

        return X[x]

    mod = Model(Modified_SEIR_deaths)

    # Initial values and bounds
    mod.set_param_hint('N', value=N, vary=False)

    mod.set_param_hint('D', value=D, vary=False)

    mod.set_param_hint('CFR_O', value=0.01, min=0.0001, max=0.1)
    mod.set_param_hint('CFR_S', value=0.1, min=0.0001, max=1.0)
    
    mod.set_param_hint('R_0_start', value=2.5, min=1.0, max=5.0)
    mod.set_param_hint('R_0_end', value=0.7, min=0.01, max=5.0)

    mod.set_param_hint('x0', value=30.0, min=0.0, max=float(Max_d))
    mod.set_param_hint('k', value=0.1, min=0.01, max=5.0)
    '''
    if R_0_kwargs:
        for arg in R_0_kwargs:
            mod.set_param_hint(arg, value=R_0_kwargs[arg])
    '''

    params = mod.make_params()
    params.add('R_0_delta', value=1.0, min=0.0, expr="R_0_start - R_0_end")  
    result = mod.fit(y_data, params, method=fit_method, x=x_data)

    CFR_O = result.params["CFR_O"].value
    CFR_S = result.params["CFR_S"].value
    R_0_result_params = {}
    for val in R_0_kwargs:
        R_0_result_params[val] = result.params[val].value

    return result, CN, y_data, N, D, Max_d, CFR_O, CFR_S, R_0_result_params

    # Comparison plot
def Modified_SEIR_fitted_plot(result, CN, y_data):
    x_ticks = pd.date_range(end=Last_Date, periods=len(y_data))

    plt.figure(figsize=(10,5))
    x_data = np.linspace(0, len(y_data), len(y_data))
    real_data, = plt.plot(x_data, y_data, 'b-', label="Reported Data")
    SIR_fit = plt.plot(x_data, result.best_fit, 'r-', label="Predicted by Present Study")
    
    plt.xlabel("Day")
    plt.ylabel("Fatalities")
    plt.legend(numpoints=1, loc=2, frameon=None)
    plt.savefig('Real Data vs SEIR-Model_NJ.png', dpi=600)    

    plt.show()


# In[11]:


result, CN, y_data, N, D, Max_d, CFR_O, CFR_S, R_0_result_params = Fit_Modified_SEIR("New_Jersey", missing_days=30, fit_method="least_squares", 
                                                                                                                 R_0_start=2.5, k=0.3, x0=170, R_0_end=0.2)
Modified_SEIR_fitted_plot(result, "New_Jersey", y_data);

plot_Modified_SEIR(*Modified_SEIR(N, D, Max_d, CFR_O, CFR_S, New_R_0, **R_0_result_params))

plot_Modified_SEIR_1(*Modified_SEIR(N, D, Max_d, CFR_O, CFR_S, New_R_0, **R_0_result_params))


# In[12]:


# Data extraction for the rest of analysis
print (Modified_SEIR(N, D, Max_d, CFR_O, CFR_S, New_R_0, **R_0_result_params))


# In[13]:


# Goodnes of fit (R2)
y=y_data
yhat = result.best_fit
SS_Residual = sum((y-yhat)**2)       
SS_Total = sum((y-np.mean(y))**2)     
r_squared = 1 - (float(SS_Residual))/SS_Total
print (r_squared) 






