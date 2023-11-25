from pathlib import Path
from PyQt5 import uic, QtWidgets
from xps_background import xps_BackgroundWindow
import numpy as np

class XPS_ApplyAllWindow(QtWidgets.QDialog):
    def __init__(self):
        super(XPS_ApplyAllWindow,self).__init__()
        self.variables = False
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
        self.variables_checkbox = self.findChild(QtWidgets.QCheckBox,"checkBox_Variables")       
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
        #window is now displayed
        self.show()
        return

    def enable_prefix(self):
        if self.save_checkbox.isChecked():
            self.prefix_lineedit.setEnabled(True)
            self.save_prefix = self.prefix_lineedit.text()
        else:
            self.prefix_lineedit.setEnabled(False)
        return

    def apply(self):
        if self.variables_checkbox.isChecked():
            self.variables = True
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




class peak_analysis:
    """Class that contains functions for peak analysis"""
    
    def apply_to_all_subwindows(self):
        """Function that applies the current settings to all subwindows"""
        #select active window
        active_window = self.mdi.activeSubWindow().widget()
        #if the active window is not none
        if active_window is not None:
            #if the active window is a graph window
            if active_window.objectName() == "PeakGraph":
                #apply settings to active window
                settings_names,settings = peak_analysis.get_settings(active_window)
                print("settings names = ",settings_names)
                #for each window in the mdi area
                for window in self.mdi.subWindowList():
                    print(window.widget().objectName(),self.mdi.subWindowList().index(window))
                    #if the window is a graph window
                    if window.widget().objectName() == "PeakGraph":
                        #if the window is not the active window
                        #apply settings to the window
                        peak_analysis.set_settings(settings_names,settings,window.widget())
                        peak_analysis.enable_GUI_functions(window.widget())
        return

    def get_settings(active_window):
        """Gets parameters from active window"""     
        settings_names = ["axis","sind","background_range","background_type","background","peak_no","inital_peak_positions","constraints","model_type","model_func","is_convoluded","model_prefix","fitting_method","model","pars","fit_results","fitted_array","peaks"]
        settings = []
        for name in settings_names:
            settings.append(getattr(active_window,name))
        return settings_names,settings
    
    def set_settings(settings_names,settings,window):
        """Sets constraints background and ranges to active window"""
        for name in settings_names:
            setting = settings[settings_names.index(name)]
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
                    print(background)
                except:
                    print("background calculation failed for ", self.mdi.subWindowList().index(window))
                peak_analysis.enable_GUI_functions(window.widget())
        return
    
    def fit_all_peaks(self):
        """Iterates over all mdi subwindow and runs the fit peaks function on each window"""
        for window in self.mdi.subWindowList():
            if window.widget().objectName() == "PeakGraph":
                try:
                    window.widget().initiate_fitting()
                except:
                    print("fitting failed for ", self.mdi.subWindowList().index(window))
                peak_analysis.enable_GUI_functions(window.widget())
        return
    
    def total_area(self):
        """Iterates over all mdi subwindow and runs the total area function on each window"""
        area = np.zeros(len(self.mdi.subWindowList()))
        energy = np.zeros(len(self.mdi.subWindowList()))
        for window in self.mdi.subWindowList():
            i = 0
            if window.widget().objectName() == "PeakGraph":
                try:   
                    area[i] = np.trapz(window.widget().data[:,1] - window.widget().background)
                    energy[i] = window.widget().meta_data["excitation_energy"]
                except:
                    print("fitting failed for ", self.mdi.subWindowList().index(window))
            i+=1
        return area, energy
    
    def save_mass_spectra(self):
        """Iterates through spectra and runs the save spectra function on each window"""
        for window in self.mdi.subWindowList():
            if window.widget().objectName() == "PeakGraph":
                window.widget().save_spectra(enmass = True,file_prefix = self.prefix,file_index = self.mdi.subWindowList().index(window))
        return




