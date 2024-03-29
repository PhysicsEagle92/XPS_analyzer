# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'XPS_main.ui'
#
# Created by: PyQt5 UI code generator 5.15.4
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.

from pathlib import Path
from PyQt5 import uic, QtWidgets
import sys
from XPS_analyzer import XPS_GraphWindow
from PeakAnalysis import XPS_ApplyAllWindow,XPS_YieldWindow, peak_analysis
from nexus_handling import get_nexus_data,get_nexus_data2
from nexusformat.nexus import nxload
from DataHandling import data_handling
from dat import dat_handling,dat_handling_I06, dat_handling_IELTS, xy_handling,csv_handling
import numpy as np
import os

class XPS_SpectrumWindow(QtWidgets.QDialog):
    def __init__(self):
        super(XPS_SpectrumWindow,self).__init__()
        self.spectrum_energies = None
        self.spectrum_id_init = None
        self.spectrum_id_final = None
        self.initUI() 

    def initUI(self):
        #loads the UI from the .ui file
        path = Path("UIWidgets")
        filename = path / "XPS_selectspectrums.ui"
        uic.loadUi(filename,self)
        self.energies_combo = self.findChild(QtWidgets.QComboBox,"comboBox_Select_Index")
        self.no_of_spectrums_combo = self.findChild(QtWidgets.QComboBox,"comboBox_No_Of_Spectrum")
        self.done_button = self.findChild(QtWidgets.QPushButton,"button_Done")
        self.energies_combo.currentIndexChanged.connect(self.load_data)
        self.no_of_spectrums_combo.currentIndexChanged.connect(self.load_data)
        #list comprehension to create a list of numbers from 0 to 49
        combonumber_list = [x for x in range(50)]
        #adds combobox list to no_of_spectrums_combo box
        self.no_of_spectrums_combo.addItems([str(x) for x in combonumber_list])
        self.done_button.clicked.connect(self.done1)
        #window is now displayed
        self.show()
        return

    def set_energy_item(self):
        self.energies_combo.clear()
        for i in range(len(self.spectrum_energies)):
            item = "ind" + str(i) + ":" + str(self.spectrum_energies[i])
            self.energies_combo.addItem(item)
        return
    
    def load_data(self):
        self.spectrum_id_init = self.energies_combo.currentIndex()
        self.spectrum_id_final = int(self.spectrum_id_init) + int(self.no_of_spectrums_combo.currentText())
        return
    
    def done1(self):
        self.close()
        return
    
class XPS_UI_Main(QtWidgets.QMainWindow):
    """Class containing all the information for the main window of program, the multiple display area (mdi)
    the functions called here load the data from source files and distribute it to windows where it needs to go"""
    
    #class constructor
    def __init__(self):
        #keeps track of the number of windows open
        self.window_count = 0
        #the file path of the file to be loaded
        self.file = None
        #the data from a vamas file, if applicable
        self.vms = None
        #the data from a dat file, if applicable
        self.dat = None
        #the data from a xy file, if applicable
        self.xy = None
        #the data from a nexus file, if applicable
        self.nxs_data = None
        #the metadata from a nexus file, if applicable
        self.nxs_meta = None
        #defines if data is loaded from a map
        self.is_mapp = False
        #the spectrum id of the first spectrum in the map
        self.spectrum_id_init = None
        #the spectrum id of the last spectrum in the map
        self.spectrum_id_final = None
        #makes the class a subclass of QMainWindow
        super(XPS_UI_Main,self).__init__()
        #loads the UI from the .ui file and assigns values to objects in that file
        self.initUI()
        
    def initUI(self):
        #loads the UI from the .ui file
        path = Path("UIWidgets")
        filename = path / "XPS_UI.ui"
        uic.loadUi(filename,self)

        #assigns values to objects, mdi area and graf button
        self.mdi = self.findChild(QtWidgets.QMdiArea, "mdiArea")
        self.graf_button = self.findChild(QtWidgets.QPushButton,"graf")
        self.grafsss = self.findChild(QtWidgets.QPushButton,"grafsss")
        #connects the button to the function that opens the graph window
        self.graf_button.clicked.connect(self.spectrum_single)
        self.grafsss.clicked.connect(self.from_map)
        #sets the menu bar in the UI MDI window
        self.set_menu_functions()
        
        #window is now displayed
        self.show()
        return

    def set_menu_functions(self):
        self.menu = self.menuBar()
        #adds file menu to toolbar
        self.file_menu = self.menu.addMenu("File")
        #adds exit to file menu
        self.Load = QtWidgets.QAction("Load Previous Spectrum", self)
        self.Load_Mass = QtWidgets.QAction("Load Previous Spectrum from directory", self)
        self.exit = QtWidgets.QAction("Exit", self)
        self.file_menu.addAction(self.Load)
        self.file_menu.addAction(self.Load_Mass)
        self.file_menu.addAction(self.exit)
        self.Load.triggered.connect(self.load_previous_spectrum)
        self.Load_Mass.triggered.connect(self.load_previous_spectrum_mass)
        self.exit.triggered.connect(self.done1)
        
        #Adds the options menu to the toolbar
        self.options_menu = self.menu.addMenu("Options")
        #adds the apply to all option to the options menu
        self.apply_all_action = QtWidgets.QAction("Apply to All...", self)
        self.options_menu.addAction(self.apply_all_action)
        self.apply_all_action.triggered.connect(self.apply_all)

        #adds the apply to all option to the options menu
        self.total_e_yield = QtWidgets.QAction("Total Electron Yield...", self)
        self.options_menu.addAction(self.total_e_yield)
        self.total_e_yield.triggered.connect(self.yield_calc)        
        return

    def done1(self):
         self.close()
         return
    
    def spectrum_single(self):
        self.is_mapp = False
        self.window_graph()
        return

    def from_map(self):
        self.is_mapp = True
        self.window_graph()
        return

    def window_graph(self):
        """Function that opens a graph window and loads the data from the file into it, nexus files can contain multiple 
        spectra so the function will open a new window for each spectrum in the file"""
        
        #loads the data from the file
        self.load_data()
        #if the file is a nexus file, it will open a new window for each spectrum in the file
        if self.file[-3:] == "nxs":
            for i, _ in enumerate(self.nxs_data):
                #adds one to window count
                self.window_count += 1
                #initilez the class containing the graph window object and its functions
                self.graph_window = XPS_GraphWindow()
                #sets the title of the window to the file name
                self.graph_window.setWindowTitle(self.file)
                #assigns the data of the ith nexus region to the graph window object
                self.graph_window.nxs = self.nxs_data[i]
                self.graph_window.meta_data = self.nxs_meta
            #initalizes the data for the region (initalize data function on the XPS_analyzer file)
                self.graph_window.initalize_data()
                #adds the graph window to the mdi area
                self.mdi.addSubWindow(self.graph_window)
                #shows the graph window
                self.graph_window.show()

        #if the file is a vamas file, it will open a new window for the file
        elif self.file[-3:] == "vms":
            self.window_count += 1
            #adds one to window count
            self.graph_window = XPS_GraphWindow()
            #sets the title of the window to the file name
            self.graph_window.setWindowTitle(self.file)
            #assigns data to graph window object
            self.graph_window.vms = self.vms
            #initalizes the data for the region (initalize data function on the XPS_analyzer file)
            self.graph_window.initalize_data()
            #adds the graph window to the mdi area
            self.mdi.addSubWindow(self.graph_window)
            #shows the graph window
            self.graph_window.show()
        #if the file is a dat file, it will open a new window for the file
        elif self.file[-3:] == "dat":
            for i, _ in enumerate(self.dat):
                self.window_count += 1
                #adds one to window count
                self.graph_window = XPS_GraphWindow()
                #sets the title of the window to the file name
                self.graph_window.setWindowTitle(self.file)
                #assigns data to graph window object
                self.graph_window.dat = self.dat[i]
                #initalizes the data for the region (initalize data function on the XPS_analyzer file)
                self.graph_window.initalize_data()
                #adds the graph window to the mdi area
                self.mdi.addSubWindow(self.graph_window)
                #shows the graph window
                self.graph_window.show()
        elif self.file[-3:] == self.file[-3:] == ".xy":
            self.window_count += 1
            #adds one to window count
            self.graph_window = XPS_GraphWindow()
            #sets the title of the window to the file name
            self.graph_window.setWindowTitle(self.file)
            #assigns data to graph window object
            self.graph_window.xy = self.xy
            #initalizes the data for the region (initalize data function on the XPS_analyzer file)
            self.graph_window.initalize_data()
            #adds the graph window to the mdi area
            self.mdi.addSubWindow(self.graph_window)
            #shows the graph window
            self.graph_window.show()       
        #if neither is possible displays this error message
        else:
            print("could not open file")
        return

    def load_data(self):
       
        """Function that opens a file dialog and loads the file name as a string into the file variable,
       if the file is a nexus file or a vamas file"""
       
       #opens a file dialog and assigns the file name to the file variable
        self.file, _ = QtWidgets.QFileDialog.getOpenFileName(self,"Open VMS File","", "Vamas Files (*.vms);;Nexus Files (*.nxs);;All Files (*)")
        # if its a nexus file calls the functions nxload and get nexus data to create a variable containing the data
        if self.file[-3:] == "nxs" and self.is_mapp == False:
            nexus_file = nxload(self.file)
            self.nxs_data, self.nxs_meta = get_nexus_data(nexus_file)
        elif self.file[-3:] == "nxs" and self.is_mapp == True:
            file = nxload(self.file)
            self.nxs_data, self.nxs_meta = get_nexus_data2(file)
            excitation_energy = self.nxs_meta[0]["excitation_energy"]           
            energies = self.nxs_meta[0]["energies"][0]
            spectrum_data = self.nxs_meta[0]["spectrum_data"]            
            #opens box to allow user to select spectrum range from map
            self.spectrum_window =  XPS_SpectrumWindow()
            self.spectrum_window.setWindowTitle("Select Spectrum")
            self.spectrum_window.spectrum_energies = excitation_energy
            self.spectrum_window.set_energy_item()
            self.spectrum_window.exec_()
            self.spectrum_id_init = self.spectrum_window.spectrum_id_init
            self.spectrum_id_final = self.spectrum_window.spectrum_id_final
            hvi = np.arange(self.spectrum_id_init,self.spectrum_id_final)
            self.nxs_data = []
            for i in range(len(hvi)):
                self.nxs_data.append({"x" : energies, "y" : spectrum_data[hvi[i]]}) 
        # if its a vamas file calls the function open vamas file to create a variable containing the data
        elif self.file[-3:] == "vms":
            self.vms = data_handling.open_vamas_file(self.file)
        #if neither is possible the function returns without any assignment
        elif self.file[-3:] == "dat":
            self.dat = dat_handling_IELTS(self.file)
        elif self.file[-3:] == ".xy":
            self.xy = xy_handling(self.file)
        elif self.file[-3:] == "csv":
            self.xy = csv_handling(self.file)
        else:
            pass
        return
    

    def apply_all(self):
        """Function that opens the apply all options window"""
        apply_all_window = XPS_ApplyAllWindow()
        apply_all_window.exec_()
        variables = apply_all_window.variables
        background = apply_all_window.background
        fit = apply_all_window.fit
        save = apply_all_window.save
        prefix = apply_all_window.save_prefix 
        if variables == True:
            peaks,contraints,ranges,backgrounds,fittings = apply_all_window.peak_settings,apply_all_window.constraint_settings,apply_all_window.datarange_settings,apply_all_window.background_settings,apply_all_window.fit_settings
            peak_analysis.apply_to_all_subwindows(self,peaks,contraints,ranges,backgrounds,fittings)
        if background == True:
            background_variables = {"range_1":self.mdi.activeSubWindow().widget().background_range["range_1"],"range_2":self.mdi.activeSubWindow().widget().background_range["range_2"],"type_1":self.mdi.activeSubWindow().widget().background_type["type_1"],"type_2":self.mdi.activeSubWindow().widget().background_type["type_2"]}
            peak_analysis.apply_background_to_all(self,background_variables)
        if fit == True:
            peak_analysis.fit_all_peaks(self)
        if save == True:
            dir_name  = QtWidgets.QFileDialog.getExistingDirectory(self,"Open Directory","/home",QtWidgets.QFileDialog.ShowDirsOnly)    
            peak_analysis.save_mass_spectra(self.mdi,dir_name,prefix)
        return

    def yield_calc(self):
        """Function that opens the yield calculator window"""
        yield_window = XPS_YieldWindow()
        #sets lenn to number of mdi windows
        
        yield_window.current_mdi = self.mdi
        yield_window.lenn = self.mdi.subWindowList().index(self.mdi.subWindowList()[-1])
        yield_window.initize_combos()
        yield_window.exec_()
        return
    
    def load_previous_spectrum(self,file):
        """Loads a spectrum from a csv file"""
        #opens a file dialog and assigns the file name to the file variable
        self.file, _ = QtWidgets.QFileDialog.getOpenFileName(self,"Open CSV File","", "CSV Files (*.csv);;All Files (*)")
        # if its a csv file calls the function open csv file to create a variable containing the data
        if self.file[-3:] == "csv":
            try:
                self.load_from_previous()
            except:
                print("Error loading file: ",self.file)
        return
    
    def load_previous_spectrum_mass(self):
        """Loads list of spectra from a directory"""
        #opens a file dialog and assigns the file name to the file variable
        directory = QtWidgets.QFileDialog.getExistingDirectory(self,"Open Directory","/home",QtWidgets.QFileDialog.ShowDirsOnly)
        #iterate over each file in directory
        for new_file in os.listdir(directory):
            full_file = directory + "/" + new_file
            self.file = full_file
            try:                
                self.load_from_previous()
                print("Loaded file: ",new_file)
            except:
                print("Error loading file: ",new_file)

        return
    
    def load_from_previous(self):
        """Loads a spectrum from a csv file"""
        energy, raw_data, fit, background, peaks, peak_data, _ = peak_analysis.load_csv_file(self.file)
        #if neither is possible the function returns without any assignment
        data = np.zeros((len(energy),2))
        data[:,0],data[:,1] = energy,raw_data
        self.spectrum_window =  XPS_GraphWindow()
        self.spectrum_window.setWindowTitle(self.file)
        self.spectrum_window.data = data
        self.spectrum_window.fitted_array = fit
        self.spectrum_window.background = background
        self.spectrum_window.show_background_button.setDisabled(False)
        self.spectrum_window.peak_no , _ = peak_analysis.find_number_of_peaks_n_params(self.file)
        self.spectrum_window.peaks = np.zeros((self.spectrum_window.peak_no,len(fit)))
        for i,key in enumerate(peaks.keys()):
            self.spectrum_window.peaks[i] = peaks[key]
        #find the indices spanning the non nan values of the data
        indicies = np.where(~np.isnan(fit))[0]
        #find the min and max of the indicies
        self.spectrum_window.sind = np.array([np.min(indicies),np.max(indicies)])
        self.spectrum_window.spanned_data = np.array([[np.nan,np.nan]]*len(self.spectrum_window.data))
        self.spectrum_window.spanned_data[self.spectrum_window.sind[0]:self.spectrum_window.sind[1]] = self.spectrum_window.data[self.spectrum_window.sind[0]:self.spectrum_window.sind[1]]           
        peak_analysis.set_peak_data(self.spectrum_window,peak_data,self.spectrum_window.peak_no)
        self.spectrum_window.show_peakposotions_button.setDisabled(False)
        self.spectrum_window.enable_fpeaks()
        #adds the graph window to the mdi area
        self.mdi.addSubWindow(self.spectrum_window)
        #shows the graph window
        self.spectrum_window.show()        
        return
    
def main_window():
    """Function that initalizes the main window object and runs the program"""
    app = QtWidgets.QApplication(sys.argv)
    XPS_MainWindow = XPS_UI_Main()
    sys.exit(app.exec_())
    return


#runs main_window
main_window()


