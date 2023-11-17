import lmfit as lm
import numpy as np
import re

import linecache
import sys

from PyQt5 import uic,QtWidgets
from PyQt5.QtCore import Qt
from pathlib import Path

from DataHandling import data_handling, XPS_PeakFindingWindow,XPS_PeakSettingWindow,XPS_DataRangeWindow
from PeakFitting import peak_fitting, XPS_ConstraintsWindow, XPS_SetPeakModel
from xps_background import xps_BackgroundWindow

#from Graph_maker import raw_data

import csv
import pandas as pd




class XPS_GraphWindow(data_handling,peak_fitting,QtWidgets.QMainWindow):
    """Class contains all functions and objects located on the graph window Gui object,
   the graph window is where the data is displayed and the user can interact with the data. The functions that do
   interact with the data are located in the DataHandling, PeakFitting and PeakAnalysis classes.
   Note the widget mplwidget is a custom widget that is used to display the graph, it is defined in the mplwidget.py file
   and loaded from the ui file directly."""
   
   #class constructor
    def __init__(self):
        #varables to contain data extracted from file
        self.vms = None
        self.nxs = None
        self.dat = None
        self.xy = None
        #the data contains a number of different x-axis types, this variable keeps track of which one is currently being used
        self.axis = "axis"
        #variable containing the data used to draw the graph
        self.data = None
        #variable containing the metadata if applicable
        self.meta_data = None
        #variable to contain the first and last indicies of the spanned data
        self.sind = np.zeros(2).astype(int)
        #Variable to contain the spanned data
        #self.spanned_data =  None
        #string containing the type of background that has been applied to the data
        self.background_range = { "range_1" : np.zeros(2).astype(int), "range_2" : np.zeros(2).astype(int)}
        self.background_type = {"type_1" : "None", "type_2" : "None",}
        #variable containing the background data
        self.background = None
        #variable containing the number of peaks being fitted
        self.peak_no = 1
        #variable containing the initial peak information
        self.inital_peak_positions = None
        #Variable containing the contraints on the peak model
        self.constraints = None
        #defines type of model being fitted
        self.model_type = "Voigt Peaks"
        #Variable that stores the function name passed into the generate model routine
        self.model_func = "VoigtModel_lmfit"
        #set is convolution is being used
        self.is_convoluded = False
        #variable to store the prefix used to name the peaks
        self.model_prefix = None
        #Sets a string to define the chosen fitting procedure
        self.fitting_method = "fit_peaks"
        #variable for the tyoe of model being used
        self.model = None
        #variable containing the peak parameters
        self.pars = None
        #variable containing the fit results in form of the lmfit fitted class
        self.fit_results = None
        #Variable for the fitted data as an array to draw the graphs
        self.fitted_array = None
        #Variable for the peak components of the fitted data
        self.peaks = None
        #function called to initialize the UI
        self.initUI_window()

    def initUI_window(self):
        """Function that loads the UI from the .ui file"""
        QtWidgets.QMainWindow.__init__(self)
        path = Path("UIWidgets")
        filename = path / "XPS_graph.ui"
        uic.loadUi(filename,self)
        
        #call function to set up the UI, i.e. each button and object contained in the window.
        self.setUIActions()

    def setUIActions(self):
        """function that assigns each button and object in the UI to a function"""

        self.actionChange_Title.triggered.connect(self.open_change_title)
        self.actionClose.triggered.connect(self.close_window)
        
        #defines combobox that sets the x-axis type for the graph, and connects it to a function that redraws the graph
        self.axis_combo = self.findChild(QtWidgets.QComboBox,"comboBox_Axis")
        self.axis_combo.currentTextChanged.connect(self.initalize_data)

        #defines a button that opens a dialog box to choose the background calculation you want to make
        self.background_button = self.findChild(QtWidgets.QPushButton,"pushButton_Background")
        self.background_button.clicked.connect(self.open_background)
        
        #defines a button that allows the user to change the data range
        self.range_button = self.findChild(QtWidgets.QPushButton,"button_DataRange")
        self.range_button.clicked.connect(self.open_data_range_window)
        
        #defines a button that opens a dialog window that allows the user to auto find peaks
        self.find_button = self.findChild(QtWidgets.QPushButton,"button_FPeaks")
        self.find_button.clicked.connect(self.open_findpeaks)

        #defines a button that opens a dialog window that allows the user to set peaks manually
        self.find_button2 = self.findChild(QtWidgets.QPushButton,"button_FPeaks_2")
        self.find_button2.clicked.connect(self.open_setpeaks)
        
        #defines a button that opens a dialog window that allows the user to set contraints
        self.constraints_button = self.findChild(QtWidgets.QPushButton,"button_Constraints")
        self.constraints_button.clicked.connect(self.open_constraints)

        #defines a button that removes the current model settings
        self.clear_model_button = self.findChild(QtWidgets.QPushButton,"pushButton_ClearModel")
        self.clear_model_button.clicked.connect(self.clear_model)
        
        
        self.crop_plot_button = self.findChild(QtWidgets.QCheckBox,"checkBox_CropPlot")
        self.crop_plot_button.clicked.connect(self.crop_plot) 
        
        #check box that when clicked will draw an x on the peak positions in the graph
        self.show_peakposotions_button = self.findChild(QtWidgets.QCheckBox,"checkBox_PeakPositions")
        self.show_peakposotions_button.setEnabled(False)
        self.show_peakposotions_button.clicked.connect(self.show_peak_positions) 

        #check box that when clicked will draw the background data on the graph
        self.show_background_button = self.findChild(QtWidgets.QCheckBox,"checkBox_Background")
        self.show_background_button.setEnabled(False)
        self.show_background_button.clicked.connect(self.show_background)

        #check box that will draw fitted curve on the graph (the sum of each peak in the model)
        self.show_fittingcurve_button = self.findChild(QtWidgets.QCheckBox,"checkBox_FittingCurve")
        self.show_fittingcurve_button.setEnabled(False)
        self.show_fittingcurve_button.clicked.connect(self.show_fitting)

        #check box that will draw the individual peaks on the graph
        self.show_fittedpeaks_button = self.findChild(QtWidgets.QCheckBox,"checkBox_FittedPeaks")
        self.show_fittedpeaks_button.setEnabled(False)
        self.show_fittedpeaks_button.clicked.connect(self.show_fitted_peaks)
        
        #button that allows fitting parameters to be carried through to next fitting
        self.keep_params_button = self.findChild(QtWidgets.QCheckBox,"checkBox_KeepParams")
        self.keep_params_button.setChecked(False)

        self.keep_constraints_button = self.findChild(QtWidgets.QCheckBox,"checkBox_KeepConstraints")
        self.keep_constraints_button.setChecked(False)

        #defines the table where the peak data is displayed
        self.peak_table = self.findChild(QtWidgets.QTableWidget,"tableWidget_Peaks")
        
        #defines combobox that allows the user to set its fitting options
        self.model_combo = self.findChild(QtWidgets.QComboBox,"comboBox_ModelOptions")
        self.model_combo.currentTextChanged.connect(self.set_model_type)

        #defines combobox that allows the user to set the fitting procedure
        self.fitting_procesudres = self.findChild(QtWidgets.QComboBox,"comboBox_FittingProcedures")
        self.fitting_procesudres.currentTextChanged.connect(self.set_fitting_method)
        
        #button used to initalise the model and accosiated parameters
        self.set_model_button = self.findChild(QtWidgets.QPushButton,"button_SelectModel")
        self.set_model_button.clicked.connect(self.set_model_and_parameters)
        self.set_model_button.setDisabled(True)


        #button used to initiate the fitting process
        self.fit_button = self.findChild(QtWidgets.QPushButton,"button_FitPeaks")
        self.fit_button.setDisabled(True)
        self.fit_button.clicked.connect(self.initiate_fitting)

        #button used to save the data to a file
        self.save_button = self.findChild(QtWidgets.QPushButton, "pushButton_SaveData")
        self.save_button.clicked.connect(self.save_spectra)
        self.save_button.setEnabled(True)
        return

    def open_change_title(self):
        """function that opens a dialog box that allows the user to change the title of the graph"""
        title_dialog = xps_titledialog()
        #executes the dialog box
        title_dialog.exec_()
        title = title_dialog.new_title
        self.setWindowTitle(title)
        return

    def close_window(self):
        self.close()
        return
    
    def initalize_data(self):
        """function that initalizes the data to be graphed from the 
        data as given from the file, this is also called when the x-axis type is changed"""
    
        #shows if vamas or nexus
        if self.vms is not None:
            #checks combobox for x-axis type
            self.axis = self.axis_combo.currentText()
            #selects the data from the vamas file
            self.data = data_handling.select_data_axes(self.vms,self.axis)
        #sets the data range to the full range of the data
        elif self.dat is not None:
            self.data = self.dat
            self.data[:,1] = self.data[:,1]*-1
        elif self.xy is not None:
            self.data = self.xy
        else:
            #defines data as the x and y data from the nexus file
            self.data = np.zeros((len(self.nxs['x']),2))
            self.data[:,0],self.data[:,1] = np.array([self.nxs['x'],self.nxs['y']])           
            #self.data[:,0] = np.flip(self.data[:,0])
            #self.datax = self.datax[self.datax[:,0].argsort()]
            #lar = len(raw_data)
            #self.data = np.zeros((lar,2))
            #self.data[:,0] = self.datax[:lar,0]
            #self.data[:,1] = self.datax[:lar,1]
            #self.data[:,1] = self.data[:,1] - raw_data/7

            
        #data_handling.normalise_data(self.data)
        #sets the data range to the full range of the data
        self.sind = [0,len(self.data[:,0])]
        self.spanned_data = np.array([[np.nan,np.nan]]*len(self.data))
        self.spanned_data[self.sind[0]:self.sind[1]] = self.data[self.sind[0]:self.sind[1]]
        self.fitted_array = np.array([np.nan]*len(self.data))
        #draws the data to the graph
        self.draw_graph()
        return
    
    def set_model_type(self):
        """function that sets the model type to the model selected in the combobox"""
        self.model_type = self.model_combo.currentText()
        return
 
    def set_fitting_method(self):
        """function that sets the fitting type to the model selected in the combobox"""
        #method_list = [,,,,"Adaptive Memory Programming for Global Optimization","Nelder-Mead","L-BFGS-B","Powell","Conjugate Gradient"]
        #method_input_list = [,',,'ampgo','nelder','lbfgs','powell','cg']
        methods_dict = {"Least Squares - Levenberg-Marquardt": 'leastsq', "Least Squares - Trust Region Reflective Method": 'least_squares',
                        "Differential Evolution" : 'differential_evolution', "Brute Force Method" : 'brute', "Basinhopping": "basinhopping", "Adaptive Memory Programming for Global Optimization": "ampgo",
                        "Nelder-Mead": "nelder", "L-BFGS-B": "lbfgs", "Powell": "powell", "Conjugate-Gradient": "cg","Newton-CG": "newton", "Cobyla" : 'cobyla', "BFGS":"bfgs",
                        "Truncated Newton": "tnc","Newton-CG trust-region": "trust-ncg","nearly exact trust-region": "trust-exact", "Newton GLTR trust-region" : "trust-krylov", 
                        "trust-region for constrained optimization" : "trust-constr" ,"Dogleg": "dogleg","Sequential Linear Squares Programming": "slsqp",
                        'Maximum likelihood via Monte-Carlo Markov Chain' : 'emcee', "Simplicial Homology Global Optimization" : "shgo", "Dual Annealing optimization" : "dual_annealing"}
        method  = self.fitting_procesudres.currentText()
        self.fitting_method = methods_dict[method]
        return
    
    def open_data_range_window(self):
        """function that allows the user to select a range of data to be graphed"""
        #calls data range function, which adds the spanning window and returns the data rangeu
        data_range_window = XPS_DataRangeWindow()
        data_range_window.data = self.data
        data_range_window.sind = [0,len(self.data[:,0])]
        data_range_window.graph_data_basic()
        data_range_window.exec_()
        self.sind = data_range_window.sind
        self.spanned_data = np.array([[np.nan,np.nan]]*len(self.data))
        self.spanned_data[self.sind[0]:self.sind[1]] = self.data[self.sind[0]:self.sind[1]]
        self.adjust_background_range()
        self.background = xps_BackgroundWindow.calculate_background(self,self.data,self.background_range['range_1'],self.background_range['range_2'],self.background_type['type_1'],self.background_type['type_2'])
        self.add_table_data_range()
        self.draw_graph()
        return  
   
    def adjust_background_range(self):
        """function that adjusts the background range to the data range"""
        if self.background_range['range_1'][0] < self.background_range['range_2'][0]:
            self.background_range['range_1'][0] = self.sind[0]
        else:
            self.background_range['range_2'][0] = self.sind[0]
        if self.background_range['range_1'][1] > self.background_range['range_2'][1]:
            self.background_range['range_1'][1] = self.sind[1]
        else:
            self.background_range['range_2'][1] = self.sind[1]
        return
    
    def open_background(self):
        """opens a dialog box that allows the user to choose the background calculation they want to make"""
        
        #initilzes background box from a UI file
        background_window = xps_BackgroundWindow()
        #imports the data into the background window
        background_window.data = self.data
        #imports the data range into the background window
        background_window.sind = np.copy(self.sind)
        background_window.background_range_1 = np.copy(self.sind)
        background_window.background_range_2 = np.copy(self.sind)
        #draws graph on background dialog box
        background_window.graph_data_basic_plus_background()
        #runs background window until it is closed
        background_window.exec_()
        #sets the background data to the data calculated in the background window
        self.background = background_window.background
        #sets indicies to those selected in the background
        self.sind = background_window.sind
        self.background_range['range_1'],self.background_range['range_2'] = background_window.background_range_1,background_window.background_range_2
        self.background_type['type_1'],self.background_type['type_2'] = background_window.background_type_1,background_window.background_type_2
        #changes data range to match background
        #self.add_table_data_range()
        #enables/disables the show background checkbox depending on if there is background data
        if self.background is None:
            self.show_background_button.setDisabled(True)
        else:
            self.show_background_button.setDisabled(False)
        #draws the data to the graph
        self.draw_graph()
        #fit peaks is only enable if there is a background or a guess
        self.enable_fpeaks()
        return
    
    def open_findpeaks(self):
        """opends a dialog box that allows the user to set variables to the scipy peak finding algorithm in order to set an inital guess for the fitting process"""
        findpeaks_window = XPS_PeakFindingWindow()
        #imports data to find peaks window
        findpeaks_window.data = self.spanned_data
        #runs find peaks window until it is closed
        findpeaks_window.exec_()
        #sets the inital peak info to the values calculated in the find peaks window
        self.inital_peak_positions = findpeaks_window.inital_peak_positions
        #sets the peak numpt to the numper of peaks found in the window
        self.peak_no = findpeaks_window.peak_no
        #enables/disables the show peak positions checkbox depending on if there is peak data
        if self.inital_peak_positions is None:
            self.show_peakposotions_button.setDisabled(True)
        else:
            self.show_peakposotions_button.setDisabled(False)
            #draws the data to the graph
            self.draw_graph()
            #adds data to the table
            self.itemise_table_from_find_peaks()
        #fit peaks is only enable if there is a background or a guess
        self.enable_fpeaks()
        return

    def open_setpeaks(self):
        """opens a dialog box that allows the user to set variables to the scipy peak finding algorithm in order to set an inital guess for the fitting process"""
        setpeaks_window = XPS_PeakSettingWindow()
        #imports data to find peaks window
        setpeaks_window.spanned_data = self.spanned_data
        setpeaks_window.data = self.data
        setpeaks_window.graph_data_basic_plus_points()
        setpeaks_window.set_inital_peak_positions_length()
        #runs find peaks window until it is closed
        setpeaks_window.exec_()
        #sets the peak numpt to the numper of peaks found in the window
        self.peak_no = setpeaks_window.peak_no
        #enables/disables the show peak positions checkbox depending on if there is peak data
        #self.setpeaks_window.inital_peak_positions = np.sort(self.setpeaks_window.inital_peak_positions)
        if setpeaks_window.inital_peak_positions is None:
            self.show_peakposotions_button.setDisabled(True)
        else:
            #sets the inital peak info to the values calculated in the find peaks window
            self.inital_peak_positions = setpeaks_window.inital_peak_positions
            self.inital_peak_positions = self.inital_peak_positions.astype(int)
            self.show_peakposotions_button.setDisabled(False)
            #draws the data to the graph
            self.draw_graph()            
            #adds data to the table
            self.itemise_table_from_find_peaks()
        #fit peaks is only enable if there is a background or a guess
        self.enable_fpeaks()
        return     
    
    def open_constraints(self):
        """opens a dialog box that allows the user to set constraints for the fitting process"""
        
        constraints_window = XPS_ConstraintsWindow()
        #imports data to find peaks window
        constraints_window.peak_no = self.peak_no
        constraints_window.inital_peak_positions = self.inital_peak_positions
        constraints_window.prefix = self.model_prefix
        constraints_window.pars = self.pars
        if self.constraints is None:
            constraints_window.constraints = self.set_constraints_variable(self.pars)
        else:
            constraints_window.constraints = self.constraints
        constraints_window.update_display_combobox()
        constraints_window.allow_add_button()
        #self.constraints_window.update_table()
        constraints_window.set_peak_number_options()
        #runs find peaks window until it is closed
        constraints_window.exec_()
        #sets the inital peak info to the values calculated in the find peaks window
        self.constraints = constraints_window.constraints
        return


    def add_table_columns_from_peak_no(self):
        """clears columns already in the table, and sets to the correct numper of columns with the peak no as their header"""
        self.peak_table.setColumnCount(0)
        if self.model_prefix is None:
            self.model_prefix = [""]*self.peak_no
        for i in range(len(self.model_prefix)):
            self.peak_table.insertColumn(i)
            self.peak_table.setHorizontalHeaderItem(i,QtWidgets.QTableWidgetItem("Peak {no}, ".format(no = i+1) + self.model_prefix[i]))
        return

    def add_table_rows_inital(self):
        """adds the inital rows to the table"""
        self.peak_table.setRowCount(8)
        self.peak_table.setVerticalHeaderItem(0,QtWidgets.QTableWidgetItem("Data Range (L)"))
        self.peak_table.setVerticalHeaderItem(1,QtWidgets.QTableWidgetItem("Data Range (R)"))
        self.peak_table.setVerticalHeaderItem(2,QtWidgets.QTableWidgetItem("Center"))
        self.peak_table.setVerticalHeaderItem(3,QtWidgets.QTableWidgetItem("Amplitude"))        
        self.peak_table.setVerticalHeaderItem(4,QtWidgets.QTableWidgetItem("Sigma"))
        self.peak_table.setVerticalHeaderItem(5,QtWidgets.QTableWidgetItem("Gamma"))
        self.peak_table.setVerticalHeaderItem(6,QtWidgets.QTableWidgetItem("fwhm"))
        self.peak_table.setVerticalHeaderItem(7,QtWidgets.QTableWidgetItem("height"))
        return
    
    def add_table_rows_from_params(self):
        self.make_param_list()
        self.peak_table.setRowCount(3 + len(self.param_list))
        table_items = ["Data Range (L)","Data Range (R)","Center"] + self.param_list
        self.peak_table.setVerticalHeaderLabels(table_items)
        return

    def make_param_list(self):
        """makes a list of the parameters in the fit results"""
        allpars = list(self.pars.keys())  
        prefixlen = len(self.model_prefix[0])
        allpars = [x[prefixlen:] for x in allpars]
        self.param_list = list(set(allpars))
        self.param_list = [cen for cen in self.param_list if "center" not in cen]
        self.param_list.sort()
        return
    
    def add_table_peak_positions(self):
        """adds the centre position of each peak to the table"""
        for i in range(len(self.model_prefix)):
            item_centre = self.cell(str(self.data[self.inital_peak_positions[i],0]))
            self.peak_table.setItem(2,i,QtWidgets.QTableWidgetItem(item_centre))
    
    def add_table_data_range(self):
        """adds the data range to the table"""
        if self.model_prefix is None and self.peak_no is not None:
            self.model_prefix = [""]*self.peak_no
        elif self.model_prefix is None and self.peak_no is None:
            self.model_prefix = [""]      
        for i in range(len(self.model_prefix)):
            item_dataL = self.cell(str(self.data[self.sind[0],0]))
            item_dataR = self.cell(str(self.data[self.sind[1]-1,0]))
            self.peak_table.setItem(0,i,QtWidgets.QTableWidgetItem(item_dataL))
            self.peak_table.setItem(1,i,QtWidgets.QTableWidgetItem(item_dataR))
        return
    
    def add_table_fitting_data(self):
        """adds the data to the table"""
        for param in self.fit_results.params:
            if "center" in param:
                i = 2
                j = int(param[len(self.model_prefix[0])-1])                
                item = self.cell(str(np.round(self.fit_results.params[param].value,5)))
                self.peak_table.setItem(i,j,QtWidgets.QTableWidgetItem(item))
            else:
                i = self.param_list.index(param[len(self.model_prefix[0]):])
                j = int(param[len(self.model_prefix[0])-1])
                item = self.cell(str(np.round(self.fit_results.params[param].value,5)))
                self.peak_table.setItem(i+3,j,QtWidgets.QTableWidgetItem(item))    
        return

    def add_table_fitting_data(self):
        """adds the data to the table"""
        for param in self.fit_results.params:
            j = int(re.findall(r'\d+',param)[0])           
            if "center" in param:
                i = 2        
                item = self.cell(str(np.round(self.fit_results.params[param].value,5)))
                self.peak_table.setItem(i,j,QtWidgets.QTableWidgetItem(item))
            else:
                i = self.param_list.index(param[len(self.model_prefix[0]):]) 
                item = self.cell(str(np.round(self.fit_results.params[param].value,5)))
                self.peak_table.setItem(i+3,j,QtWidgets.QTableWidgetItem(item))    
        return

    def itemise_table_from_find_peaks(self):
        self.add_table_columns_from_peak_no()
        self.add_table_peak_positions()
        self.add_table_data_range()
        return
    
    def cell(self,var):
            item = QtWidgets.QTableWidgetItem()
            item.setText(var)
            item.setFlags(item.flags() ^ Qt.ItemIsEditable)
            return item

    def show_peak_positions(self):
        """graphs data if checkbox is clicked"""
        self.draw_graph()
        return
   
    def show_background(self):
        """graphs data if checkbox is clicked"""
        self.draw_graph()
        return
    
    def show_fitting(self):
        """graphs data if checkbox is clicked"""
        self.draw_graph()
        return

    def show_fitted_peaks(self):
        """graphs data if checkbox is clicked"""
        self.draw_graph()
        return
    
    def crop_plot(self):
        """graphs data if checkbox is clicked"""
        self.draw_graph()
        return
  
    def draw_graph(self):
        if self.crop_plot_button.isChecked():
            self.graph_data_cropped()
        else:
            self.graph_data()
        return
    
    def clear_model(self):
        """clears the model"""
        if self.model == None:
            return
        else:
            self.model = None
            self.fit_results = None
            self.pars = None
            self.model_prefix = None
            self.constraints = None
            self.model_type = "Voigt Peaks"
            self.model_func = "VoigtModel_lmfit"
            self.peak_no = 1
            self.peak_table.clear()
            self.add_table_columns_from_peak_no()
            self.add_table_rows_inital()            
            self.add_table_data_range()
            self.find_button.setEnabled(True)
            self.find_button2.setEnabled(True)
        return
    
    def set_model(self):
        if self.model_type == "Voigt Peaks":
            self.model_func = "VoigtModel_lmfit"
        elif self.model_type == "Lorentzian Peaks":
            self.model_func = "LorentzianModel_lmfit"
        elif self.model_type == "Doniach-Sunjic Peaks":
            self.model_func = "DoniachModel_lmfit"
        elif self.model_type == "Convolved Doniach-Sunjic and Guassian":
             self.model_func = "Convolved_DoniachGuassian_lmfit"
             self.is_convoluded = True
        elif self.model_type == "reversed Doniach-Sunjic Peak":
            self.model_func = "reversed_DoniachModel_lmfit"
        elif self.model_type == "reversed Convolved Doniach-Sunjic and Guassian":
            self.model_func = "reveresed_Convolved_DoniachGuassian_lmfit"
            self.is_convoluded = True
        elif self.model_type == "Set Individual Peaks":
            self.set_peak_model = XPS_SetPeakModel(self.peak_no,self.model_prefix)
            self.set_peak_model.exec_()
            self.model_func = self.set_peak_model.set_model
            self.model,self.model_prefix = peak_fitting.set_individual_peaks(self,self.peak_no,self.model_func)
            return
        else:
            print("unable to select model")
            return

        set_model = getattr(peak_fitting,self.model_func)
        self.model,self.pars,self.model_prefix = set_model(self,self.peak_no,self.spanned_data[self.inital_peak_positions,0])
        return
    
    def set_model_and_parameters(self):
        self.set_model()
        #checks if user wants to keep paramerters generated in the last fitting
        if self.keep_params_button.isChecked() == True and self.inital_peak_positions is not None:
            self.pars = self.keep_parameters_fixed(self.pars,self.constraints)
        self.enable_fpeaks()
        self.add_table_columns_from_peak_no()
        self.add_table_data_range()
        self.add_table_peak_positions()
        self.find_button.setEnabled(False)
        self.find_button2.setEnabled(False)
        return
    
    def get_prefix_from_constraints(self):
        """Works out peak prefic from the parameters variable"""
        keys = list(self.constraints.keys())
        prefix_check = []
        prefix_ind = keys[0].find("0")
        prefix_check = list(set([key[:prefix_ind+1] for key in keys]))
        return prefix_check
    
    def initiate_fitting(self):
        """runs the fitting algorithm and set a number of buttons to be enabled or disabled depending on the outcome"""
        #creates a generic name 'set model' which is then used to call the correct model function from a string
        x,y = self.spanned_data[:,0],self.spanned_data[:,1]-self.background
        try:
            if self.constraints is not None:
                print("from conts", self.get_prefix_from_constraints())
                print("from pres",self.model_prefix)
                if all(self.get_prefix_from_constraints()) == all(self.model_prefix):
                    print("Adding Constraints")
                    self.set_constraint_parameters()
                elif self.keep_constraints_button.isChecked() == True:
                    print("Adding with transfer")
                    self.convert_constraints_to_new_model()
                    print("done")
                    self.set_constraint_parameters()
                else:
                    print("No added constraints")
                    pass
            else:
                print("No added constraints")
                pass
        except:
            print("\nUnable to add Constraints\n")
            self.PrintException()
        #runs the fitting algorithm
        try:
            self.fit_results = peak_fitting.fit_peaks(self,x,y,self.pars,self.model,self.fitting_method)
            print("\nFitting completed!\n")
            data_handling.set_graph_variables(self)
            self.show_fittingcurve_button.setEnabled(True)
            self.show_fittedpeaks_button.setEnabled(True)
            self.display_fitting()
            self.draw_graph()
        except Exception as e: 
            print("\nFitting failed!\n")
            self.PrintException()
        #if fitting has been completed successfully the will enable the checkbox to show the fitting curve and fitted peaks, and the save peak data button
        return

    def display_fitting(self):
        """prints out fit report """
        #self.display_fit_report()
        self.add_table_rows_from_params()
        self.add_table_fitting_data()
        return

    def enable_fpeaks(self):
        """checks if there is a background or an inital peak positions and enables the fit peaks button if there is"""
        if self.inital_peak_positions is not None:
            self.set_model_button.setDisabled(False)
        else:
            self.set_model_button.setDisabled(True)

        if self.background is not None and self.inital_peak_positions is not None and self.model is not None and self.pars is not None:
            self.fit_button.setDisabled(False)
        else:
            self.fit_button.setDisabled(True)
        
        return
    
    
    def save_report(self,file_name):
        """saves infromation from the peak model into a JSON file, as provided by lmfit """        
        with open(file_name, 'w', newline='') as file:
            print("Saving Report")
            writer = csv.writer(file)
            if self.fit_results is not None:
                for keys in self.fit_results.params:
                    writer.writerow([keys,self.fit_results.params[keys].value])                    
                writer.writerow(["Chi Squared",self.fit_results.chisqr])
                writer.writerow(["Reduced Chi Squared",self.fit_results.redchi])
                writer.writerow(["AIC",self.fit_results.aic])
                writer.writerow(["BIC",self.fit_results.bic])
                writer.writerow(["#############DATA#############"])
                headers = ["Energy","Raw_Data","Fit","Background"]
                for p in range(len(self.peaks)):
                    headers.append("Peaks{}".format(p+1))
                writer.writerow(headers)
                for i in range(len(self.data)):
                    data_line = [self.data[i,0],self.data[i,1],self.fitted_array[i],self.background[i]]
                    for p in range(len(self.peaks)):
                        data_line.append(self.peaks[p][i])
                    writer.writerow(data_line)                    
            else:
                print("writing lines")
                writer.writerow(["#############DATA#############"])
                headers = ["Energy","Raw_Data"]
                writer.writerow(headers)
                print(len(self.data))
                for i in range(len(self.data)):
                    data_line = [self.data[i,0],self.data[i,1]]
                    writer.writerow(data_line)  
                    print(data_line)
        return
    
    def save_spectra(self,enmass = False,file_prefix = "Spectra_",file_index = "00"):
        if enmass == False:
            file_name, _  = QtWidgets.QFileDialog.getSaveFileName(self,"Save csv File","", "csv file (*.csv)")
            if file_name != "":
                self.save_report(file_name)
            else:
                pass
        else:
            file_name = file_prefix + str(file_index) + ".csv"
            if file_name != "":
                self.save_report(file_name)
            else:
                pass
        return
    
    def PrintException(self):
        exc_type, exc_obj, tb = sys.exc_info()
        f = tb.tb_frame
        lineno = tb.tb_lineno
        filename = f.f_code.co_filename
        linecache.checkcache(filename)
        line = linecache.getline(filename, lineno, f.f_globals)
        print('EXCEPTION IN ({}, LINE {} "{}"): {}'.format(filename, lineno, line.strip(), exc_obj))
        return
    
    
class xps_titledialog(QtWidgets.QDialog):
    def __init__(self):
        super(xps_titledialog,self).__init__()
        self.new_title = None
        self.initUI()
        
    def initUI(self):
        #loads the UI from the .ui file
        path = Path("UIWidgets")
        filename = path / "XPS_title.ui"
        uic.loadUi(filename,self)
        self.edit_title = self.findChild(QtWidgets.QLineEdit,"lineEdit_Title")
        self.done_button = self.findChild(QtWidgets.QPushButton,"button_Done")
        self.done_button.clicked.connect(self.get_new_title)
        return
    
    def get_new_title(self):
        self.new_title = self.edit_title.text()
        self.close()
        return
