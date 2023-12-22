from pathlib import Path
from PyQt5 import uic, QtWidgets
from xps_background import xps_BackgroundWindow
from DataHandling import data_handling
import numpy as np
import copy
import csv

class XPS_ApplyAllWindow(QtWidgets.QDialog):
    def __init__(self):
        super(XPS_ApplyAllWindow,self).__init__()
        self.variables = False
        self.peak_settings = False
        self.constraint_settings = False
        self.datarange_settings = False
        self.background_settings = False
        self.fit_settings = False                                                                            
        self.background = False
        self.fit = False
        self.save = False
        self.save_prefix = "spectra_"
        self.initUI() 

    def initUI(self):
        #loads the UI from the .ui file
        path = Path("UIWidgets")
        filename = path / "XPS_applyall.ui"
        uic.loadUi(filename,self)
        #sets apply button
        self.apply_button = self.findChild(QtWidgets.QPushButton,"button_Apply")
        self.apply_button.clicked.connect(self.apply)
        #defines the checkboxs of variables to apply
        self.allvariables_checkbox = self.findChild(QtWidgets.QCheckBox,"checkBox_AllVariables")
        self.allvariables_checkbox.clicked.connect(self.enable_variables)
        self.peaksettings_checkbox = self.findChild(QtWidgets.QCheckBox,"checkBox_Variable1_PeakSettings")
        self.contraints_checkbox = self.findChild(QtWidgets.QCheckBox,"checkBox_Variable2_Contraints")
        self.datarange_checkbox = self.findChild(QtWidgets.QCheckBox,"checkBox_Variable3_DataRange")
        self.background_setting_checkbox = self.findChild(QtWidgets.QCheckBox,"checkBox_Variable4_Background")
        self.fittings_checkbox = self.findChild(QtWidgets.QCheckBox,"checkBox_Variable5_Fitting")
        self.peaksettings_checkbox.setEnabled(False)
        self.contraints_checkbox.setEnabled(False)
        self.datarange_checkbox.setEnabled(False)
        self.background_setting_checkbox.setEnabled(False)
        self.fittings_checkbox.setEnabled(False)
            
        self.background_checkbox = self.findChild(QtWidgets.QCheckBox,"checkBox_BackgroundCalc")
        self.fit_checkbox = self.findChild(QtWidgets.QCheckBox,"checkBox_Fit")
        self.save_checkbox = self.findChild(QtWidgets.QCheckBox,"checkBox_Save")
        self.save_checkbox.clicked.connect(self.enable_prefix)
        #done button
        self.done_button = self.findChild(QtWidgets.QPushButton,"button_Done")
        self.done_button.clicked.connect(self.done1)
        #fins prefix lineedit
        self.prefix_lineedit = self.findChild(QtWidgets.QLineEdit,"lineEdit_SetPrefix")
        self.prefix_lineedit.setText("Spectra_")
        self.prefix_lineedit.setDisabled(True)
        self.prefix_lineedit.textChanged.connect(self.enable_prefix)
        #window is now displayed
        self.show()
        return

    def enable_variables(self):
        if self.allvariables_checkbox.isChecked():
            self.peaksettings_checkbox.setEnabled(True)
            self.contraints_checkbox.setEnabled(True)
            self.datarange_checkbox.setEnabled(True)
            self.background_setting_checkbox.setEnabled(True)
            self.fittings_checkbox.setEnabled(True)
        else:
            self.peaksettings_checkbox.setEnabled(False)
            self.contraints_checkbox.setEnabled(False)
            self.datarange_checkbox.setEnabled(False)
            self.background_setting_checkbox.setEnabled(False)
            self.fittings_checkbox.setEnabled(False)

        return
    
    def enable_prefix(self):
        if self.save_checkbox.isChecked():
            self.prefix_lineedit.setEnabled(True)
            self.save_prefix = self.prefix_lineedit.text()
        else:
            self.prefix_lineedit.setEnabled(False)
        return

    def apply(self):
        if self.allvariables_checkbox.isChecked():
            self.variables = True
            if self.peaksettings_checkbox.isChecked():
                self.peak_settings = True
            if self.contraints_checkbox.isChecked():
                self.constraint_settings = True
            if self.datarange_checkbox.isChecked():
                self.datarange_settings = True
            if self.background_checkbox.isChecked():
                self.background_settings = True
            if self.fittings_checkbox.isChecked():
                self.fit_settings = True
            
        if self.background_checkbox.isChecked():
            self.background = True
        if self.fit_checkbox.isChecked():
            self.fit = True
        if self.save_checkbox.isChecked():
            self.save = True
        self.close()
        return

    
    def done1(self):
        self.close()
        return


class XPS_YieldWindow(QtWidgets.QDialog):
    def __init__(self):
        super(XPS_YieldWindow,self).__init__()
        self.current_mdi = None
        self.lenn = None
        self.initUI() 
        

    def initUI(self):
        #loads the UI from the .ui file
        path = Path("UIWidgets")
        filename = path / "XPS_totalyield.ui"
        uic.loadUi(filename,self)
        self.start_combobox = self.findChild(QtWidgets.QComboBox,"comboBox_Start")
        self.end_combobox = self.findChild(QtWidgets.QComboBox,"comboBox_End")     
        #sets apply button
        self.apply_button = self.findChild(QtWidgets.QPushButton,"button_Apply")
        self.apply_button.clicked.connect(self.apply)
        #save button
        self.save_button = self.findChild(QtWidgets.QPushButton,"pushButton_SaveData")
        self.save_button.clicked.connect(self.save_as_csv)
        #done button
        self.done_button = self.findChild(QtWidgets.QPushButton,"button_Done")
        self.done_button.clicked.connect(self.done1)
        #window is now displayed
        self.show()
        return

    def initize_combos(self):
        """Function that initializes the comboboxesby settng them to be numbers from index 0 to the length of the spectrum"""
        numbers = [x for x in range(self.lenn)]
        self.start_combobox.addItems([str(x) for x in numbers])
        self.end_combobox.addItems([str(x) for x in numbers])
        #set value of start combo to first value
        self.start_combobox.setCurrentIndex(0)
        #set value of end combo to last value
        self.end_combobox.setCurrentIndex(self.lenn-1)
        return

    def apply(self):
        self.start = int(self.start_combobox.currentText())
        self.finish = int(self.end_combobox.currentText())
        self.eyield,self.energy = peak_analysis.total_area(self.current_mdi,self.start,self.finish)
        #draws graph of yield vs energy
        self.graph_data()
        return
    
    def graph_data(self):
        """Function to plot chosen XPS plots data set"""
        #initalizes variables X & Y for plotting data
        X,Y = self.energy,self.eyield
        self.MplWidget.canvas.axes.clear()
        self.MplWidget.canvas.axes.plot(X,Y)
        self.MplWidget.canvas.draw()
        return
    
    def save_as_csv(self):
        """Function that saves the data as a .csv file"""
        #gets the file name from the user
        filename = QtWidgets.QFileDialog.getSaveFileName(self,"Save File","C:\\","CSV (*.csv)")
        #if the filename is not empty
        if filename[0] != "":
            #if the filename does not end in .csv
            if filename[0][-4:] != ".csv":
                #add .csv to the end of the filename
                filename = filename[0] + ".csv"
            else:
                filename = filename[0]
            #writes the data to the file
            with open(filename,"w") as f:
                f.write("Energy,Yield\n")
                for i in range(len(self.energy)):
                    f.write(str(self.energy[i]) + "," + str(self.eyield[i]) + "\n")
        return
    
    def done1(self):
        self.close()
        return

class peak_analysis:
    """Class that contains functions for peak analysis"""
    
    def apply_to_all_subwindows(self,peaks,contraints,ranges,backgrounds,fittings):
        """Function that applies the current settings to all subwindows"""
        #select active window
        active_window = self.mdi.activeSubWindow().widget()
        #if the active window is not none
        if active_window is not None:
            #if the active window is a graph window
            if active_window.objectName() == "PeakGraph":
                #apply settings to active window
                settings_names,settings = peak_analysis.get_settings(active_window)
                #for each window in the mdi area
                for window in self.mdi.subWindowList():
                    #if the window is a graph window
                    if window.widget().objectName() == "PeakGraph":
                        #if the window is not the active window
                        #apply settings to the window
                        peak_analysis.set_settings(settings_names,settings,window.widget(),peaks,contraints,ranges,backgrounds,fittings)
                        peak_analysis.enable_GUI_functions(window.widget())
        return

    def get_settings(active_window):
        """Gets parameters from active window"""     
        settings_names = {"peak_settings" : ["peak_no","inital_peak_positions","model_type","is_convoluded","model_func","model_prefix","fitting_method","model","pars"],"constraints":["constraints"],"data_range":["sind"],"background" : ["background_range","background_type"],"fittings":["fit_results","fitted_array","peaks"]}        
        settings = {"peak_settings" : [],"constraints":[],"data_range":[],"background":[],"fittings":[]}
        for settings_key in settings_names.keys():
            for name in settings_names[settings_key]:  
                settings[settings_key].append(getattr(active_window,name))
        return settings_names,settings
    
    def set_settings(settings_names,settings,window,peaks,contraints,ranges,backgrounds,fittings):
        """Sets constraints background and ranges to active window"""
        setting_keys = []
        if peaks == True:
            setting_keys.append("peak_settings")
        if contraints == True:
            setting_keys.append("constraints")
        if ranges == True:
            setting_keys.append("data_range")
        if backgrounds == True:
            setting_keys.append("background")
        if fittings == True:
            setting_keys.append("fittings")
        for key in setting_keys:
            for name in settings_names[key]:
                i = settings_names[key].index(name)
                setting = copy.deepcopy(settings[key][i])
                setattr(window,name,setting)
        return
    
    def enable_GUI_functions(window):
        """Enables GUI functions for window"""
        #updates if the fitting buttons are allowed
        window.enable_fpeaks()
         #updates if the background buttons are allowed
        if window.background is None:
            window.show_background_button.setDisabled(True)
        else:
            window.show_background_button.setEnabled(True)
        #updates if the peak positions buttons are allowed
        if window.inital_peak_positions is None:
            window.show_peakposotions_button.setDisabled(True)
        else:
            window.show_peakposotions_button.setDisabled(False)
        #updates the tabel if a mpdel is present
        if window.model is None:
            pass
        else:
            window.add_table_columns_from_peak_no()
            window.add_table_data_range()
            window.add_table_peak_positions()
            window.find_button.setEnabled(False)
            window.find_button2.setEnabled(False)
        #updates if a fitting is able to be displayed
        if window.fit_results is None:
            window.show_fittingcurve_button.setDisabled(True)
            window.show_fittedpeaks_button.setDisabled(True)           
        else:
            window.show_fittingcurve_button.setEnabled(True)
            window.show_fittedpeaks_button.setEnabled(True)
            window.display_fitting()
    
        return
 
    def apply_background_to_all(self,background_variables):
        """iterates over all subwindows runniong the background calculation"""
        for window in self.mdi.subWindowList():
            if window.widget().objectName() == "PeakGraph":
                background = xps_BackgroundWindow.calculate_background(window.widget(),window.widget().data,background_variables['range_1'],background_variables['range_2'],background_variables['type_1'],background_variables['type_2'])              
                try:
                    background = xps_BackgroundWindow.calculate_background(window.widget(),window.widget().data,background_variables['range_1'],background_variables['range_2'],background_variables['type_1'],background_variables['type_2'])
                    window.widget().background = background
                except:
                    print("background calculation failed for ", self.mdi.subWindowList().index(window))
                peak_analysis.enable_GUI_functions(window.widget())
        return
    
    def fit_all_peaks(self):
        """Iterates over all mdi subwindow and runs the fit peaks function on each window"""
        #prints the index of the last window in mdi subwindow list
        for window in self.mdi.subWindowList():
            print("Now fitting window", self.mdi.subWindowList().index(window) + 1, " of ",len(self.mdi.subWindowList()))
            if window.widget().objectName() == "PeakGraph":
                try:
                    window.widget().initiate_fitting()
                except:
                    print("fitting failed for ", self.mdi.subWindowList().index(window))
        return
    
    def total_area(mdi,start,finish):
        """Iterates over all mdi subwindow and runs the total area function on each window"""
        window_list = peak_analysis.create_list(start,finish)
        area = np.zeros(len(window_list))
        energy = np.zeros(len(window_list))
        for window_id in window_list:
            window = mdi.subWindowList()[window_id]
            if window.widget().objectName() == "PeakGraph":    
                area[window_id] = np.trapz(window.widget().data[:,1])
                energy[window_id] = window.widget().meta_data[0]["excitation_energy"][window_id]
        return area, energy
    
    def create_list(a,b):
        if a>b:
            return [x for x in range(b,a+1)]
        else:       
            return [x for x in range(a,b+1)]    
        
    def save_mass_spectra(mdi,dir_name,prefix):
        """Iterates through spectra and runs the save spectra function on each window"""
        
        for window in mdi.subWindowList():
            if window.widget().objectName() == "PeakGraph":
                try:
                    window.widget().save_spectra(enmass = True,dir_name = dir_name,file_prefix = prefix,file_index = mdi.subWindowList().index(window))
                except:
                    print("save spectra failed for ", mdi.subWindowList().index(window))
        return
    
    
    def load_csv_file(filename):
        """Reads CSV files from XPS analyzer and returns data in NumPy arrays."""
        with open(filename, 'r') as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=',')

            # Initialize data arrays
            energy, raw_data, fit, background = np.array([]), np.array([]), np.array([]), np.array([])

            # Initialize dictionaries for header and peak/fit data
            header_data, peak_data, fit_data, peaks = {}, {}, {}, {}
        
            # Variables
            ite,peak_params = peak_analysis.find_number_of_peaks_n_params(filename)
            fit_params = 4
            find_data_line = False
            for i in range(ite):
                peaks["peak_"+str(i+1)] = np.array([])       

            
            for row in csv_reader:
                # Header lines
                if not find_data_line and row[0][0] != "#":
                    header_data[row[0]] = row[1]
                elif not find_data_line and row[0][0] == "#":
                    find_data_line = True
               
                
                # Data lines
                elif find_data_line:
                    if row[0] == "Energy":
                        pass
                    else:
                        # Add data to arrays
                        energy = np.append(energy, float(row[0]))
                        raw_data = np.append(raw_data, float(row[1]))
                        fit = np.append(fit, float(row[2]))
                        background = np.append(background, float(row[3]))
                    
                        # Add peak data to the 'peaks' dictionary
                        for i in range(ite):                        
                            peaks["peak_"+str(i+1)] = np.append(peaks["peak_"+str(i+1)],float(row[i+4]))   
                    
        # Split header data into 'peak_data' and 'fit_data'
        
            param_counter = 0
            for key, value in header_data.items():            
                if param_counter < peak_params*ite:
                    peak_data[key] = value
                else:
                    fit_data[key] = value
                param_counter +=1
        return energy, raw_data, fit, background, peaks, peak_data, fit_data


    
    def find_number_of_peaks_n_params(filename):
        """Finds the number of peaks in the CSV file."""
        with open(filename, 'r') as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=',')
            no_of_rows = 0
            for row in csv_reader:
                no_of_rows += 1
                #if statement that reads true if the row starts with a letter
                if row[0][0].isalpha():
                    pass
                elif row[0] == "Energy":
                    pass
                elif row[0][0] == "#":
                    pass
                else:
                    peak_no = len(row) - 4
                    no_params = int((no_of_rows -3 - 4)/peak_no)
                    return peak_no, no_params
            
    def set_peak_data(spectrum_window,peak_data,peak_no):
        """unpacks the peak_data and converts it to the spectrum window"""
        model_types_dict = {"V" : "Voigt Peaks" , "L" : "Lorentzian Peaks", "D" : "Doniach-Sunjic Peaks", "rD" : "reversed Doniach-Sunjic Peak", 
                            "C" : "Convolved Doniach-Sunjic and Guassian", "rC" :"reversed Convolved Doniach-Sunjic and Guassian"}
        #get model type from peak data
        key = next(iter(peak_data))
        spectrum_window.model_type = model_types_dict[key[0]]
        #get all keys that end with "center"
        center_keys = [key for key in peak_data.keys() if key.endswith('center')]
        peak_positions = np.array([],dtype = 'int')
        for key in center_keys:           
            peak_index = np.abs(spectrum_window.data[:,0] - float(peak_data[key])).argmin()
            peak_positions = np.append(peak_positions,peak_index)
        spectrum_window.inital_peak_positions = peak_positions
        spectrum_window.set_model_and_parameters()
        for key in peak_data.keys():
            if key.endswith('center'):
                pass
            else:
                spectrum_window.pars[key].set(value = float(peak_data[key]))
        spectrum_window.itemise_table_from_find_peaks()
        spectrum_window.add_table_data_range()
        spectrum_window.show_fittingcurve_button.setEnabled(True)
        spectrum_window.show_fittedpeaks_button.setEnabled(True)
        spectrum_window.make_param_list()
        spectrum_window.add_table_fitting_data_from_pars()
        spectrum_window.draw_graph()
        #run set model function in spectrum window
        
        return


