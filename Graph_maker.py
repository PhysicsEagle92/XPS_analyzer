import numpy as np
import matplotlib.pyplot as plt
import csv
from matplotlib.animation import FuncAnimation
from xps_background import background_calculation
import lmfit as lm
import os
#import data from csv file
def import_csv(filename):
    """function to read csv files from XPS analyzer and return data in numpy arrays"""
    #opens csv file
    with open(filename) as csv_file:
        #reads csv file
        data = csv.reader(csv_file, delimiter=',')
        #creates empty arrays to store data
        #electron kinetic energy array (x-axis)
        energy = np.array([])
        #raw data variable (y-axis)
        raw_data = np.array([])
        #fit data variable
        fit = np.array([])
        #background data variable
        background = np.array([])
        #number of stored parameters for each peak
        peak_params = 6
        #number of stored parameters related to the fit
        fit_params = 4
        #dictionary to store header data
        header_data = {}
        #dictionary to store peak data
        peak_data = {}
        #dictionary to store fit data
        fit_data = {}
        #dictionary to store data for each peak
        peaks = {}
        #variable to determine when data starts
        find_data_line = 0
        #variable to count lines
        line_count = 0
        #iterates through each line in the csv file
        for row in data:
            #checks if the line is a header line
            if row[0][0] != "#" and find_data_line == 0:
                #adds header data to header_data dictionary
                header_data[row[0]] = row[1]
                line_count += 1
            #csv deliminates header data and graph data by row of hashs (mostly)
            elif row[0][0] == "#" and find_data_line == 0:
                find_data_line = 1
                line_count += 1
                #works out number of parameters for fit and peak data so the header data can be split into the correct dictionaries, which is equal to the number of peaks
                ite = 6#int((line_count-fit_params - 1)/peak_params) + 1
                #print(ite)
                #iterates over the number of header lines and adds the data to the correct dictionary
                for i in range(ite):
                    #also defines the number of peaks
                    peaks["peak_"+str(i+1)] = np.array([])
                    count = 0
                    for keys in header_data:
                        #adds the peak data to the peak_data dictionary
                        if count < peak_params*ite:
                            peak_data[keys] = header_data[keys]
                            count+=1
                        else:
                            #adds the fit data to the fit_data dictionary
                            fit_data[keys] = header_data[keys]
            #checks if we have reaced the headers for the data columns
            elif row[0] == "Energy" and find_data_line == 1:
                pass
            #checks if we have reached the end of the data
            elif find_data_line == 1:
                #adds data to the correct arrays
                energy = np.append(energy,float(row[0]))
                raw_data = np.append(raw_data,float(row[1]))
                fit = np.append(fit,float(row[2]))
                background = np.append(background,float(row[3]))
                #iterates over the number of peaks and adds the data to the correct arrays
                for i in range(ite):
                    peaks["peak_"+str(i+1)] = np.append(peaks["peak_"+str(i+1)],float(row[i+4]))                
                line_count += 1
            
    return energy,raw_data,fit,background,peaks,peak_data,fit_data

def animate_points(energy,raw_data,fit,peaks, range):
    """funtion to graph the peak data, raw data and fit """
    fig, axs = plt.subplots(1,1)
    def animate(i):
        axs.clear()
        axs.scatter(energy[range],raw_data[range],marker = '.',color = 'blue')
        axs.plot(energy[range],fit[range],color = 'red')
        axs.scatter(energy[i],raw_data[i],marker = 'X',color = 'green')
        axs.set_ylabel('Intensity (arb.uni.)')
        for keys in peaks:
            axs.fill(energy[range],peaks[keys][range],alpha = 0.5)
        axs.legend(['raw data','fit'])
        axs.set_xlabel('Energy (eV)')
        return
    ani = FuncAnimation(fig, animate, frames=len(y), interval=50, repeat=False)
    plt.show()
    return

def open_two_column_csv(file):
    """function to open a csv file with two columns and one header line and return the data as two numpy arrays"""
    with open(file) as csv_file:
        data = csv.reader(csv_file, delimiter=',')
        x = np.array([])
        y = np.array([])
        line_count = 0
        for row in data:
            if line_count == 0:
                line_count += 1
            else:
                x = np.append(x,float(row[0]))
                y = np.append(y,float(row[1]))
                line_count += 1
    return x,y
    


def graph_with_residels(energy,raw_data,fit,background,peaks, range):
    """funtion to graph the peak data, raw data and fit """
    fig, axs = plt.subplots(2,1,sharex = True,gridspec_kw={'height_ratios': [3, 1]})
    fig.subplots_adjust(hspace=0)
    axs[0].scatter(energy[range],raw_data[range]-background[range],marker = '.',color = 'blue')
    axs[0].plot(energy[range],fit[range]-background[range],color = 'red')
    axs[0].set_ylabel('Intensity (a.u.)')    
    print(peaks.keys())
    i = 0
    for keys in peaks:
        if i < 3:
            axs[0].fill(energy[range],peaks[keys][range]-background[range],facecolor = "black",alpha = 0.6)
        else:
            axs[0].fill(energy[range],peaks[keys][range]-background[range],facecolor = "grey",alpha = 0.6)
        i += 1
    axs[0].legend(['raw data','fit'])
    axs[1].set_ylabel('percent dif (%)')
    axs[0].set_xlabel('Energy (eV)')
    axs[1].scatter(energy,(1-(fit/raw_data))*100,marker = '.',color = 'blue')
    plt.show()
    return
    
def graph_without_residels(energy,raw_data,fit,peaks, range):
    """funtion to graph the peak data, raw data and fit """
    fig, axs = plt.subplots(1,1)
    axs.scatter(energy[range],raw_data[range],marker = '.',color = 'blue')
    axs.plot(energy[range[0:22]],fit[range[0:22]],color = 'red')
    axs.set_ylabel('Intensity (a.u.)',fontsize = 16)
    for keys in peaks:
        axs.fill(energy[range],peaks[keys][range],alpha = 0.5)
    axs.legend(['raw data','fit'])
    axs.set_xlabel('Energy (eV)',fontsize = 16)
    plt.show()
    return

def load_all_peaks(file_list):
    peaks = []
    fits = []
    for name in file_list:
        data = import_csv(name)
        fits.append(data[2])
        peaks.append(data[3])
    return data[0],data[1],fits,peaks

def stitch_peaks_graph(file_list):
    """funtion to graph the peak data, raw data and fit """
    energy,raw_data,fit,peaks = load_all_peaks(file_list)
    background = linear_background(energy,raw_data)
    fig, axs = plt.subplots(1,1)
    axs.scatter(energy,(raw_data-background)*-1,marker = '.',color = 'blue')
    #axs.plot(energy,background)
    axs.set_ylabel('Intensity (a.u.)',fontsize = 16)
    axs.set_xlabel('Energy (eV)',fontsize = 16)
    for i in range(len(peaks)):
        nan_range = find_nan_range(energy,fit[i])
        axs.plot(energy[nan_range],(fit[i][nan_range]-background[nan_range])*-1,color = 'red')
        for keys in peaks[i]:
            axs.fill(energy[nan_range],(peaks[i][keys][nan_range]-background[nan_range])*-1,alpha = 0.5)
    plt.show()
    return

def linear_background(energy,raw_data):
	"""Calculates a background by subtracting a linear line set between two points"""
	p1,p2 = np.array([energy[0],raw_data[0]]),np.array([energy[-1],raw_data[-1]])
	m = (p2[1]-p1[1])/(p2[0]-p1[0])
	b = p1[1]-m*p1[0]
	background = np.array([m*x+b for x in energy])
	return background
    
def find_nan_range(energy,peak):
    """function to find the range of data that is not nan (the nans are from changing the range in the XPS analyzer"""
    nan_range = np.array([]).astype(int)   
    for i in range(len(energy)):
        if np.isnan(peak[i]) == False:
            nan_range = np.append(nan_range,i)
    return nan_range


def write_csv(filename,energy,raw_data):
    """function to write data to a csv file"""
    with open(filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        writer.writerow(['Energy','Raw Data'])
        for i in range(len(energy)):
            writer.writerow([energy[i],raw_data[i]])
    return

            
def calcculate_area_under_lorenztian(amplitude,gamma):
    """function to calculate the area under a lorenztian curve"""
    area = float(amplitude)*float(gamma)*np.pi
    return area


def charge_transfer_v_photon_energy(filelist,photon_energy,chl):
    """function to calculate the charge transfer from the area under the peak"""
    ct = np.zeros(len(filelist))
    i = 0
    for fil in filelist:
        print(fil)
        data = import_csv(fil)
        peak_data = data[4]
        auger_area = calcculate_area_under_lorenztian(peak_data["L0amplitude"],peak_data["L0sigma"]) + calcculate_area_under_lorenztian(peak_data["L1amplitude"],peak_data["L1sigma"]) + calcculate_area_under_lorenztian(peak_data["L2amplitude"],peak_data["L2sigma"])
        res_area = calcculate_area_under_lorenztian(peak_data["L3amplitude"],peak_data["L3sigma"]) + calcculate_area_under_lorenztian(peak_data["L4amplitude"],peak_data["L4sigma"]) + calcculate_area_under_lorenztian(peak_data["L5amplitude"],peak_data["L5sigma"])
        
        ct[i] = (res_area/auger_area)*chl
        i+=1
    #exp_model = lm.models.ExponentialModel()
    #pars = exp_model.guess(ct,x=photon_energy)
    #pars["decay"].set(value = 1,vary=True,min=0.01,max=100)
    #pars["amplitude"].set(value = 10,vary=True,min=0.01,max = 10000)
    #print(pars["decay"])
    #print(pars["amplitude"])
    #out = exp_model.fit(ct,pars,x=photon_energy-245,method = "differential_evolution")
    #print(out.fit_report(min_correl=0.25))
    #x = np.arange(0,7,0.01)/100

    #plt.plot(x+245,out.eval(x=x))
   #plt.plot(photon_energy,out.best_fit)
    plt.scatter(photon_energy,ct)
    plt.xlabel('Photon Energy (eV)')
    plt.ylabel('Charge Transfer (fs)')
    plt.show()
    return

           
def get_peak_area(peaks,ranges):
    """function to integrate peak area"""
    #calculates integral from peak data
    area = np.array([]).astype(float) 
    for key in peaks:
        area = np.append(area,np.trapz(peaks[key][ranges]))
    return area

def total_area(data):
    return np.trapz(data)

def draw_lorentz(x,x0,A,g):
    return (A/np.pi)*((0.5*g)/((x-x0)**2 + (0.5*g)**2))

def cal_comp():
    E_calc_gas = np.array([400.53,400.81,401.11,401.42,401.72])
    A_meas_gas = np.array([188.452,184.590,96.317,35.717,8.973])
    E_calc_C60 = np.array([400.48,400.78,401.09,401.40,401.70])
    A_meas_C60 = np.array([188.452,184.590,96.317,35.717,8.973])
    x = np.arange(energy[0],energy[-1],0.0001)
    y = np.zeros((len(E_calc_gas),len(x))) 
    num_points = len(x)
    num_peaks = len(E_calc_gas)
    peak_width = 0.1368

    background = linear_background(energy,raw_data)

    for i in range(num_peaks):
        y[i] = draw_lorentz(x,E_calc_C60[i],A_meas_gas[i],peak_width)
        #plt.fill(x,y[i],color = 'blue',alpha = 0.5)
    y = np.sum(y,axis=0)*scale
    line1, = plt.plot(x,y,color = 'orange')
    line2, = plt.plot(energy-0.06+scalew,raw_data-background,color = 'blue')
    plt.legend(['Calculation','Experiment'])

    plt.show()
    return


def animate_graphs_from_directory(directory):
    """function to animate graphs from a directory"""
    files = os.listdir(directory)
    files.sort()
    file_data = []
    for fil in files:
        if fil.endswith(".csv"):
            data = import_csv(directory+fil)
            file_data.append(data)
    def animate(i):
        plt.clf()
        plt.scatter(file_data[i][0],file_data[i][1]-file_data[i][3])
        #plt.plot(file_data[i][0],file_data[i][2]-file_data[i][3],color = 'red')
        for peak in range(3):            
            plt.fill(file_data[i][0],file_data[i][4]["peak_"+str(peak+1)]-file_data[i][3],color = 'black',alpha = 0.7)
        for peak in range(3):
            plt.fill(file_data[i][0],file_data[i][4]["peak_"+str(peak+4)]-file_data[i][3],color = 'grey',alpha = 0.7)
        fit = sums_all_non_NAN_values_in_a_dictionary(file_data[i][4])
        plt.plot(file_data[i][0],fit-file_data[i][3]*6,color = 'red')
        plt.xlabel('Energy (eV)')
        plt.ylabel('Intensity (arb. units)')
        plt.title(files[i])
    ani = FuncAnimation(plt.gcf(),animate,frames=len(files),interval=300)
    #save gif animation
    ani.save("C:\\Users\\ppxcf1\\Desktop\\Data\\X@C60\\I09-beamtime-November\\Processed Data\\Ar@C60_XPS.gif", writer='imagemagick', fps=1)
    plt.show()
    return

def calcualte_peak_area(directory):
    """function to calculate peak area from a directory"""
    files = os.listdir(directory)
    files.sort()
    area_ratio = []
    for fil in files:
        peak_data_Aug = []
        peak_data_Res = []
        if fil.endswith(".csv"):
            data = import_csv(directory+fil)
            background = data[3]
            for peak in range(3):
                peak_data_Aug.append(data[4]["peak_"+str(peak+1)])
            for peak in range(3):
                peak_data_Res.append(data[4]["peak_"+str(peak+4)])
            peak_area_Aug = peak_area_calc(peak_data_Aug,background)
            peak_area_Res = peak_area_calc(peak_data_Res,background)            
            area_ratio.append(6*np.sum(peak_area_Res)/np.sum(peak_area_Aug))
    #area_ratio.remove(area_ratio[12])
    #area_ratio.remove(area_ratio[22])
    #area_ratio.remove(area_ratio[32])
    #area_ratio.remove(area_ratio[41])
    #area_ratio.remove(area_ratio[41])
    #area_ratio.remove(area_ratio[41])
    #area_ratio.remove(area_ratio[41])
    #area_ratio.remove(area_ratio[41])
    x = np.arange(1,len(area_ratio)+1,1)
    fit = fit_with_exponetial(x,area_ratio)
    print(fit.fit_report())
    plt.plot(x,fit.best_fit)
    plt.scatter(x,area_ratio)
    plt.show()
    return

def fit_with_exponetial(x,y):
    exp_model = lm.models.ExponentialModel()
    pars = exp_model.guess(y,x=x)
    pars["decay"].set(value = 1,vary=True,min=0.01,max=100)
    pars["amplitude"].set(value = 10,vary=True,min=0.01,max = 10000)
    print(pars["decay"])
    print(pars["amplitude"])
    out = exp_model.fit(y,pars,x=x,method = "differential_evolution")
    print(out.fit_report(min_correl=0.25))
    #x = np.arange(0,7,0.01)/100

    #plt.plot(x+245,out.eval(x=x))
   #plt.plot(photon_energy,out.best_fit)
    return out

def peak_area_calc(peak_data,background):
    """function to calculate peak area"""
    peak_area = np.zeros(len(peak_data))
    for i in range(len(peak_data)):
        peak_data[i] = peak_data[i][~np.isnan(peak_data[i])] - background[~np.isnan(background)]
        peak_area[i] = np.trapz(peak_data[i])
    return peak_area

def sums_all_non_NAN_values_in_a_dictionary(dict):
    """function to sum all non NAN values in a dictionary"""
    sum = 0
    for key in dict:
        #extracts non NAN values from array
        #dict[key] = dict[key][~np.isnan(dict[key])]
        sum += dict[key]
    return sum

##################C60_XPS####################
csv_file_N2aC60 = 'C:\\Users\\ppxcf1\\Desktop\\Data\\X@C60\\I09-beamtime-January\\Processed_Data\\237114_Nitrogen_XAS peak_fitted.csv'
csv_file_N2aC60_best = 'C:\\Users\\ppxcf1\\Desktop\\Data\\I09-beamtime-January\\Processed_Data\\237114_Nitrogen_XAS peak_fitted_best_fit.csv'
csv_file_N2gas = 'C:\\Users\\ppxcf1\\Desktop\\Data\\X@C60\\I09-beamtime-January\\Processed_Data\\N2_gas_phase_peaks.csv'
csv_file_Ar2P = 'C:\\Users\\ppxcf1\\Desktop\\Data\\I09-beamtime-January\\Processed_Data\\236669_Ar2P_ArMulti_Voigt_fitted peaks.csv'
csv_file_C1S_ArMulti = 'C:\\Users\\ppxcf1\\Desktop\\Data\\I09-beamtime-January\\Processed_Data\\236667_C1s_ArMulti_Voigt_fitted_peaks.csv'
csv_file_C1S_KrMulti = 'C:\\Users\\ppxcf1\\Desktop\\Data\\I09-beamtime-January\\Processed_Data\\236708_C1s_KrMulti_Voigt_fitted_peaks.csv'
csv_file_Kr3p = 'C:\\Users\\ppxcf1\\Desktop\\Data\\I09-beamtime-January\\Processed_Data\\236720_Kr3p_Voigt_fitted_peaks.csv'
csv_file_C1S_N2Multi = 'C:\\Users\\ppxcf1\\Desktop\\Data\\I09-beamtime-January\\Processed_Data\\237099_C1s_N2Multi_Voigt_fitted_peaks.csv'
#csv_file_C1S_COMulti = 'C:\\Users\\ppxcf1\\Desktop\\Data\\I09-beamtime-January\\Processed_Data\\237147_C1s_COMulti_Voigt_fitted_peaks.csv'
csv_file_Kr3d_KrMulti = 'C:\\Users\\ppxcf1\\Desktop\\Data\\I09-beamtime-January\\Processed_Data\\236735_Kr3d_KrMulti_Voigt_fitted_peaks.csv'
csv_file_Ar2p_ArMono = 'C:\\Users\\ppxcf1\\Desktop\\Data\\I09-beamtime-January\\Processed_Data\\237160_Ar2p_ArMono_Doniach_fitted_peaks.csv'
csv_file_N1s_N2Multi = 'C:\\Users\\ppxcf1\\Desktop\\Data\\I09-beamtime-January\\Processed_Data\\237103_N1s_N2Multi_Voigt_fitted_peaks.csv'
csv_file_C1s_ArMono = 'C:\\Users\\ppxcf1\\Desktop\\Data\\I09-beamtime-January\\Processed_Data\\237157_C1s_ArMono_DoniachnVoigt_fitted_peaks.csv'
#sv_file_C1s_ArMono2 = 'C:\\Users\\ppxcf1\\Desktop\\Data\\I09-beamtime-January\\Processed_Data\\237157_C1s_ArMono_DoniachnVoigt_fitted_peaks.csv'
#csv_file_Ar2p_ArMono = 'C:\\Users\\ppxcf1\\Desktop\\Data\\I09-beamtime-January\\Processed_Data\\2236986_Overview_no_fitting.csv'
csv_file_N2aC60_broad_XAS = 'C:\\Users\\ppxcf1\\Desktop\\Data\\I09-beamtime-January\\Processed_Data\\237113_Nitrogen_XAS_broad_view_voigt_fitted.csv'
csv_file_N2aC60_fullbroad_XAS = 'C:\\Users\\ppxcf1\\Desktop\\Data\\I09-beamtime-January\\Processed_Data\\237113_Nitrogen_XAS_fullbroad_view_voigt_fitted.csv'
csv_file_AraC60_XAS = 'C:\\Users\\ppxcf1\\Desktop\\Data\\I09-beamtime-January\\Processed_Data\\236692_ArXAS_ArMulti_Lorenenztian_fitted_peak.csv'

#ResPes
csv_file_Ar_ResPes_AugerPeak_offres_for_sub = 'C:\\Users\\ppxcf1\\Desktop\\Data\\I09-beamtime-January\\Processed_Data\\236699_Ar_ResPes_AugerPeak_offres_for_sub.csv'
csv_file_AraC60_ResPes_Auger_peak1 = 'C:\\Users\\ppxcf1\\Desktop\\Data\\I09-beamtime-January\\Processed_Data\\236696_Ar_ResPes_AugerPeak1_209p49eV.csv'
csv_file_AraC60_ResPes_offres = "C:\\Users\\ppxcf1\\Desktop\\Data\\I09-beamtime-January\\Processed_Data\\236696_Ar_ResPes_AugerPeak_offres.csv"
csv_file_AraC60_ResPes_offresm3 = "C:\\Users\\ppxcf1\\Desktop\\Data\\I09-beamtime-January\\Processed_Data\\Ar_ResPes_SpectrumNo20_OffRes_m3.csv"
csv_file_AraC60_ResPes_offresm2 = "C:\\Users\\ppxcf1\\Desktop\\Data\\I09-beamtime-January\\Processed_Data\\Ar_ResPes_SpectrumNo21_OffRes_m2.csv"
csv_file_AraC60_ResPes_offresm1 = "C:\\Users\\ppxcf1\\Desktop\\Data\\I09-beamtime-January\\Processed_Data\\Ar_ResPes_SpectrumNo22_OffRes_m1.csv"
csv_file_AraC60_ResPes_onres = "C:\\Users\\ppxcf1\\Desktop\\Data\\I09-beamtime-January\\Processed_Data\\Ar_ResPes_SpectrumNo23_OnRes3.csv"
csv_file_AraC60_ResPes_offresp1 = "C:\\Users\\ppxcf1\\Desktop\\Data\\I09-beamtime-January\\Processed_Data\\Ar_ResPes_SpectrumNo24_OffResp1.csv"
csv_file_AraC60_ResPes_offresp2 = "C:\\Users\\ppxcf1\\Desktop\\Data\\I09-beamtime-January\\Processed_Data\\Ar_ResPes_SpectrumNo25_OffResp2.csv"
csv_file_AraC60_ResPes_offresp3 = "C:\\Users\\ppxcf1\\Desktop\\Data\\I09-beamtime-January\\Processed_Data\\Ar_ResPes_SpectrumNo26_OffResp3.csv"
csv_Ar_ResPes_Mono_SpectrumNo14_OnRes = "C:\\Users\\ppxcf1\\Desktop\\Data\\I09-beamtime-January\\Processed_Data\\Ar_ResPes_Mono_SpectrumNo14_OnRes.csv"
csv_Ar_ResPes_Mono_SpectrumNo13_OffResm1 = "C:\\Users\\ppxcf1\\Desktop\\Data\\I09-beamtime-January\\Processed_Data\\Ar_ResPes_Mono_SpectrumNo13_OffResm1.csv"
csv_co_XAS = "C:\\Users\\ppxcf1\\Desktop\\Data\\I09-beamtime-January\\Processed_Data\\237143_CO_XAS_vibrations.csv"
csv_co_XAS_237143 = "C:\\Users\\ppxcf1\\Desktop\\Data\\X@C60\\I09-beamtime-January\\Processed_Data\\237143_CO_XAS_vibrations.csv"

csv_file_ARaC60_ResPes_Multilayer_Nov2023_Dir = "C:\\Users\\ppxcf1\\Desktop\\Data\\X@C60\\I09-beamtime-November\\Processed Data\\254771_MultiLayer_Voigt\\"
csv_file_ARaC60_ResPes_Multilayer_Nov2023_Dir2 = "C:\\Users\\ppxcf1\\Desktop\\Data\\X@C60\\I09-beamtime-November\\Processed Data\\I09-254771_Multilayer_Voigt_Improved\\"
csv_file_ARaC60_ResPes_Multilayer_Nov2023_Dir3 = "C:\\Users\\ppxcf1\\Desktop\\Data\\X@C60\\I09-beamtime-November\\Processed Data\\I09-254771_Multilayer_Voigt_AnG_fixed\\"
csv_file_ARaC60_ResPes_TEY = "C:\\Users\\ppxcf1\\Desktop\\Data\\X@C60\\I09-beamtime-November\\Processed Data\\254771_XAS_TotalElectronYield.csv"
csv_file_ARaC60_ResPes_OnRes = "C:\\Users\\ppxcf1\\Desktop\\Data\\X@C60\\I09-beamtime-November\\Processed Data\\254771_OnResSpectrum.csv"
diff = [0.251]#eV

#Nanoparticles
csv_file_CoNp_For_Quantity_Comparison = 'C:\\Users\\ppxcf1\\Desktop\\Data\\I06-beamtime-April\\CoNp_For_Quantity_Comparison.csv'
csv_file_FeNp_For_Quantity_Comparison_badfitting = 'C:\\Users\\ppxcf1\\Desktop\\Data\\I06-beamtime-April\\FeNp_For_Quantity_Comparison_badfitting.csv'
csv_file_checkshape= 'C:\\Users\\ppxcf1\\Desktop\\Data\\I06-beamtime-April\\XAS_check_shape.csv'
csv_file_CoNp_PeakFitting = 'C:\\Users\\ppxcf1\\Desktop\\Data\\I06-beamtime-April\\CoNp_PeakFitting.csv'
csv_file_CoNp_PeakFitting = 'C:\\Users\\ppxcf1\\Desktop\\Data\\I06-beamtime-April\\CoNp_PeakFitting_gamma.csv'
csv_file_CoNp_PeakFitting2 = 'C:\\Users\\ppxcf1\\Desktop\\Data\\I06-beamtime-April\\CoNp_PeakFitting_gamma2.csv'


    

##################IELTSC60#####################
csv_file_KrC60IELTSV = "C:\\Users\\ppxcf1\\Desktop\\Data\\Unisoku\\IELTS_C60\\Processed_Data\\C60-Pb-April_00325_m0p3476_voight_fitted.csv"
csv_file_KrC60IELTSL = "C:\\Users\\ppxcf1\\Desktop\\Data\\Unisoku\\IELTS_C60\\Processed_Data\\C60-Pb-April_00325_m0p3476_lorenzian_fitted.csv"
csv_file_KrC60IELTS_Hg1 = "C:\\Users\\ppxcf1\\Desktop\\Data\\Unisoku\\IELTS_C60\\Processed_Data\\C60-Pb-April_00086_Hg1_lorenzian_fitted.csv"
csv_file_KrC60IELTS_Hg2 = "C:\\Users\\ppxcf1\\Desktop\\Data\\Unisoku\\IELTS_C60\\Processed_Data\\C60-Pb-April_00086_Hg2_lorenzian_fitted.csv"
csv_file_KrC60IELTS_Hg3Hg4 = "C:\\Users\\ppxcf1\\Desktop\\Data\\Unisoku\\IELTS_C60\\Processed_Data\\C60-Pb-April_00086_Hg3Hg4_lorenzian_fitted.csv"
csv_file_KrC60IELTS_Hg5 = "C:\\Users\\ppxcf1\\Desktop\\Data\\Unisoku\\IELTS_C60\\Processed_Data\\C60-Pb-April_00086_Hg5_lorenzian_fitted.csv"
csv_file_KrC60IELTS_Hg6 = "C:\\Users\\ppxcf1\\Desktop\\Data\\Unisoku\\IELTS_C60\\Processed_Data\\C60-Pb-April_00086_Hg6_lorenzian_fitted.csv"
csv_file_KrC60IELTS_Hg7 = "C:\\Users\\ppxcf1\\Desktop\\Data\\Unisoku\\IELTS_C60\\Processed_Data\\C60-Pb-April_00086_Hg7_lorenzian_fitted.csv"
csv_file_KrC60IELTS_Hg8 = "C:\\Users\\ppxcf1\\Desktop\\Data\\Unisoku\\IELTS_C60\\Processed_Data\\C60-Pb-April_00086_Hg8_lorenzian_fitted.csv"

#################MnTeXPS########################
csv_file_Mn2p_grazing_set1_unreal = "C:\\Users\\ppxcf1\\Desktop\\Data\\MnTe\\MnTe_XPS\\ProcessedData\\Grazing_70_21_Mn2p_set1_unreal.csv"
csv_file_Mn2p_grazing_set2_unreal = "C:\\Users\\ppxcf1\\Desktop\\Data\\MnTe\\MnTe_XPS\\ProcessedData\\Grazing_70_21_Mn2p_set2_unreal.csv"
#print("spacing = ",287.43296-287.16391)
#file_list = [csv_file_KrC60IELTS_Hg1,csv_file_KrC60IELTS_Hg2,csv_file_KrC60IELTS_Hg3Hg4,csv_file_KrC60IELTS_Hg5,csv_file_KrC60IELTS_Hg6,csv_file_KrC60IELTS_Hg7,csv_file_KrC60IELTS_Hg8]

#file_list = [csv_file_AraC60_ResPes_offresm3,csv_file_AraC60_ResPes_offresm2,csv_file_AraC60_ResPes_offresm1,csv_file_AraC60_ResPes_onres,csv_file_AraC60_ResPes_offresp1,csv_file_AraC60_ResPes_offresp2,csv_file_AraC60_ResPes_offresp3]
#file_list = [csv_Ar_ResPes_Mono_SpectrumNo13_OffResm1,csv_Ar_ResPes_Mono_SpectrumNo14_OnRes]
#photon_energy = np.arange(0,7)/100)+245

#stuff = import_csv(csv_file_N2gas)
#peak_area = total_area(stuff[2])
#energy = stuff[0]
#raw_data = stuff[1]
#charge_transfer_v_photon_energy(file_list,photon_energy,6)
#ranger = find_nan_range(energy,stuff[3]['peak_1'])
#data_with_back = stuff[2][ranger] - background_calculation._calculate_shirley_background_full_range(stuff[2][ranger])
#area = np.trapz(data_with_back)
#print(area)
#graph_with_residels(energy,stuff[1],stuff[2],stuff[3],ranger)
#graph_without_residels(energy,stuff[1],stuff[2],stuff[3],ranger)
#stitch_peaks_graph(file_list)

#Fearea = total_area(data_with_back)
#print(Fearea)

#stuff = import_csv(csv_file_ARaC60_ResPes_OnRes)
#energy = stuff[0]
#ranger = find_nan_range(energy,stuff[4]['peak_1'])
#graph_with_residels(energy,stuff[1],stuff[2],stuff[3],stuff[4],ranger)

#p1 = stuff[3]['peak_1'][ranger] - background_calculation._calculate_shirley_background_full_range(stuff[1][ranger])
#p2 = stuff[3]['peak_2'][ranger] - background_calculation._calculate_shirley_background_full_range(stuff[1][ranger])
#p3 = stuff[3]['peak_3'][ranger] - background_calculation._calculate_shirley_background_full_range(stuff[1][ranger])
#p4 = stuff[3]['peak_4'][ranger] - background_calculation._calculate_shirley_background_full_range(stuff[1][ranger])
#Coarea = total_area(data_with_back)

#a1 = total_area(p1) 
#a2 = total_area(p2) 
#a3 = total_area(p3) 
#a4 = total_area(p4) 
#print("ratio of Co to CoO is: ", a2/(a1+a4+a3))
#stuff[1][ranger]  = (stuff[1][ranger] - np.min(stuff[1][ranger]))
#stuff[1][ranger] = stuff[1][ranger]/np.max(stuff[1][ranger])
#write_csv('C:\\Users\\ppxcf1\\Desktop\\Data\\I09-beamtime-January\\Processed_Data\\N2aC60_normalized.csv',energy[ranger],stuff[1][ranger])
scale = 1#922.88/1707.98


#x,y =  open_two_column_csv(csv_file_ARaC60_ResPes_TEY)
#plt.plot(x,y)
#plt.show()
#animate_graphs_from_directory(csv_file_ARaC60_ResPes_Multilayer_Nov2023_Dir2)
calcualte_peak_area(csv_file_ARaC60_ResPes_Multilayer_Nov2023_Dir3)

#make list of files in a directory
#filelist = []
#for file in os.listdir(csv_file_ARaC60_ResPes_Multilayer_Nov2023_Dir2):
#    if file.endswith(".csv"):
#        filelist.append(os.path.join(csv_file_ARaC60_ResPes_Multilayer_Nov2023_Dir2, file))
        
#photon_energy = np.arange(1,len(filelist)+1,1)
#charge_transfer_v_photon_energy(filelist,photon_energy,6)