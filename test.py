import numpy as np
from matplotlib import rc
import matplotlib.pyplot as plt
from scipy.constants import h, c
from morse import Morse, FAC
from nexusformat.nexus import nxload
from scipy.signal import fftconvolve
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from scipy.integrate import simps
from scipy.special import wofz
import lmfit as lm
from nexus_handling import get_nexus_data_spectra
from nexusformat.nexus import nxload
from xps_background import background_calculation

rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size': 14})
rc('text', usetex=True)

def morse_plot():
    fig, ax = plt.subplots()
    X.plot_V(ax, color='white')

    X.draw_Elines(range(X.vmax), ax)
    X.draw_Elines(X.get_vmax(), ax, linestyles='--', linewidths=2)
    #X.plot_psi([2], ax, scaling=2, color=COLOUR1)
    X.label_levels([0,1,2,3,4], ax)

    ax.set_xlabel(r'$r\;/\mathrm{\\A}$')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    fig.patch.set_alpha(0.0)
    fig.patch.set_facecolor('grey')
    ax.patch.set_facecolor('grey')
    ax.patch.set_alpha(0.0)
    plt.savefig('morse-psi.png')
    plt.show()
    
#morse_plot()

#fig, ax = plt.subplots()
#X.plot_psi_map([0,1,2,3,4], ax,h/FAC, scaling=2,color=COLOUR1)
#ax.patch.set_alpha(0.0)    
#fig.patch.set_facecolor('grey')
#ax.patch.set_facecolor('grey')
#ax.patch.set_alpha(0.0)  
#plt.show()
def linear_background(energy,raw_data):
	"""Calculates a background by subtracting a linear line set between two points"""
	p1,p2 = np.array([energy[0],raw_data[0]]),np.array([energy[-1],raw_data[-1]])
	m = (p2[1]-p1[1])/(p2[0]-p1[0])
	b = p1[1]-m*p1[0]
	background = np.array([m*x+b for x in energy])
	return background

def VoigtModel_lmfit2(self,peak_no):
    """Function that creates an object array, each representing one peak"""
    #initializes an empty array, before populating each element with a VoightModel object
    #Each peak is labeled in the fit report by the prefix v{i}
    peaks = np.empty(peak_no,dtype = "object")
    model_prefix = ["V{}".format(i) for i in range(peak_no)]
    for i in range(peak_no):
        peaks[i] = lm.models.VoigtModel(prefix = model_prefix[i])
    return peaks,model_prefix

def VoigtModel_lmfit(peak_no,center_pos):
    """Function that creates an object array, each representing one peak"""
    #initializes an empty array, before populating each element with a VoightModel object
    #Each peak is labeled in the fit report by the prefix v{i}
    peaks = np.empty(peak_no,dtype = "object")
    model_prefix = ["V{}".format(i) for i in range(peak_no)]
    for i in range(peak_no):
        peaks[i] = lm.models.VoigtModel(prefix = model_prefix[i])
    model = np.sum(peaks)
    pars = model.make_params()
    for i in range(len(peaks)):
        pars[model_prefix[i] + "amplitude"].set(value = 10,min = 1)
        pars[model_prefix[i]+"center"].set(center_pos[i],min = center_pos[i]-0.1,max = center_pos[i]+0.1)
        pars[model_prefix[i]+"sigma"].set(value = 0.1,min = 0.01,max = 0.5)
        pars[model_prefix[i]+"gamma"].set(value = 0.1,min = 0.01,max = 0.5,expr = None,vary = True)
    return model,pars,model_prefix

def test_model():
    file_name = "C:\\Users\\ppxcf1\\Desktop\\Data\\X@C60\\I09-beamtime-January\\si31574\\si31574-1\\i09-236735.nxs"
    nxdata = nxload(file_name)
    data_region_list, metadata_region_list = get_nexus_data_spectra(nxdata)
    x = data_region_list[0]['x']
    y = data_region_list[0]['y']
    background = linear_background(x,y)
    y = y - background
    center_pos = np.array([92.34,93.61])

    model,pars,model_prefix = VoigtModel_lmfit(2,center_pos)

    #peak1 = lm.models.VoigtModel(prefix = "V0")
    #peak2 = lm.models.VoigtModel(prefix = "V1")
    #model = peak1 + peak2
    #pars = model.make_params()
    #pars['V0center'].set(center_pos[0],min = 92.2,max = 92.4)
    #pars['V1center'].set(center_pos[1],min = 93.5,max = 93.7)
    pars['V1amplitude'].set(expr = 'V0Ampltiude*7')
    result = model.fit(y,pars,x=x)
    print(pars['V0center'].value)
    print(pars['V1center'].value)
    #print(result.fit_report())
    plt.figure()
    plt.plot(x,result.best_fit)
    plt.scatter(x,y)
    plt.show()
    return


def sin_func(x,a):
    return np.sin(x*a)

def poly(x,a,b,c):
    return a*x**2 + b*x + c


def line(x,m,c):
    return m*x + c

def plot_function(x,y):
    plt.figure()
    plt.plot(x,y)
    return


from wrapt_timeout_decorator import timeout
import time


# Define a decorator with a timeout of 5 seconds

def timeout_function(x):
    # Your long-running code here
    time.sleep(x)
    return

@timeout(5,use_signals=False)
def test_timeout(x):
    try:
        timeout_function(x)
    except:
        print("Timeout")
        return 1
    return 0

#test1 = 0
test1 = test_timeout(10)
print(test1)

test2 = test_timeout(1)
print(test2)
    



    

