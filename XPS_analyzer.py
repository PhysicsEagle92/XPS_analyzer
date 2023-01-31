import lmfit as lm
import numpy as np

from PyQt5 import uic,QtWidgets
from PyQt5.QtCore import Qt
from pathlib import Path

from DataHandling import data_handling, XPS_PeakFindingWindow,XPS_PeakSettingWindow,XPS_DataRangeWindow
from PeakFitting import peak_fitting
from xps_background import xps_BackgroundWindow


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
        #the data contains a number of different x-axis types, this variable keeps track of which one is currently being used
        self.axis = "axis"
        #variable containing the data used to draw the graph
        self.data = None
        #variable containing the metadata if applicable
        self.meta_data = None
        #variable to contain the first and last indicies of the spanned data
        self.sind = np.zeros(2).astype(int)
        #string containing the type of background that has been applied to the data
        self.background_type = "None"
        #variable containing the background data
        self.background = None
        #variable to contain the first and last indicies of the spanned data of background
        self.gind = np.zeros(2).astype(int)
        #variable containing the number of peaks being fitted
        self.peak_no = 1
        #variable containing the initial peak information
        self.inital_peak_info = None
        #variable for the tyoe of model being used
        self.model = None
        #variable containing the peak parameters
        self.pars = None
        #variable containing the fit results
        self.fit_results = None
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
        
        #defines the table where the peak data is displayed
        self.peak_table = self.findChild(QtWidgets.QTableWidget,"tableWidget_Peaks")
        
        #button used to initiate the fitting process
        self.fit_button = self.findChild(QtWidgets.QPushButton,"button_FitPeaks")
        self.fit_button.setDisabled(True)
        self.fit_button.clicked.connect(self.initiate_fitting)

        #button used to save the data to a file
        self.save_button = self.findChild(QtWidgets.QPushButton, "pushButton_SaveData")
        self.save_button.clicked.connect(self.save_report)
        self.save_button.setEnabled(False)
        return

    def open_change_title(self):
        """function that opens a dialog box that allows the user to change the title of the graph"""
        self.title_dialog = xps_titledialog()
        #executes the dialog box
        self.title_dialog.exec_()
        title = self.title_dialog.new_title
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
        else:
            #defines data as the x and y data from the nexus file
            self.data = np.zeros((len(self.nxs['x']),2))
            self.data[:,0],self.data[:,1] = np.array([self.nxs['x'],self.nxs['y']])
        #sets the data range to the full range of the data
        self.sind = [0,len(self.data[:,0])]
        self.gind = [0,len(self.data[:,0])]
        #draws the data to the graph
        self.graph_data()
        return

    def open_background(self):
        """opens a dialog box that allows the user to choose the background calculation they want to make"""
        
        #initilzes background box from a UI file
        self.background_window = xps_BackgroundWindow()
        #imports the data into the background window
        self.background_window.data = self.data
        #runs background window until it is closed
        self.background_window.exec_()
        #sets the background type to the type selected in the background window
        self.background_type = self.background_window.background_type
        #sets the background data to the data calculated in the background window
        self.background = self.background_window.background
        #enables/disables the show background checkbox depending on if there is background data
        if self.background is None:
            self.show_background_button.setDisabled(True)
        else:
            self.show_background_button.setDisabled(False)
        #draws the data to the graph
        self.graph_data()
        #fit peaks is only enable if there is a background or a guess
        self.enable_fpeaks()
        return

    def open_findpeaks(self):
        """opends a dialog box that allows the user to set variables to the scipy peak finding algorithm in order to set an inital guess for the fitting process"""
        self.findpeaks_window = XPS_PeakFindingWindow()
        #imports data to find peaks window
        self.findpeaks_window.data = self.data[self.sind[0]:self.sind[1],:]
        #runs find peaks window until it is closed
        self.findpeaks_window.exec_()
        #sets the inital peak info to the values calculated in the find peaks window
        self.inital_peak_info = self.findpeaks_window.inital_peak_info
        #sets the peak numpt to the numper of peaks found in the window
        self.peak_no = self.findpeaks_window.peak_no
        #enables/disables the show peak positions checkbox depending on if there is peak data
        if self.inital_peak_info is None:
            self.show_peakposotions_button.setDisabled(True)
        else:
            self.show_peakposotions_button.setDisabled(False)
            #draws the data to the graph
            self.graph_data()
            #corrects for spanning in data indicies
            for i in range(self.peak_no):
                self.inital_peak_info[0][i] = self.inital_peak_info[0][i] + self.sind[0]
            #adds data to the table
            self.itemise_table_from_find_peaks()
        #fit peaks is only enable if there is a background or a guess
        print(self.inital_peak_info)
        self.enable_fpeaks()
        return

    def open_setpeaks(self):
        """opends a dialog box that allows the user to set variables to the scipy peak finding algorithm in order to set an inital guess for the fitting process"""
        self.setpeaks_window = XPS_PeakSettingWindow()
        #imports data to find peaks window
        self.setpeaks_window.data = self.data[self.sind[0]:self.sind[1],:]
        self.setpeaks_window.graph_data_basic_plus_points()
        #runs find peaks window until it is closed
        self.setpeaks_window.exec_()
        #sets the peak numpt to the numper of peaks found in the window
        self.peak_no = self.setpeaks_window.peak_no
        #enables/disables the show peak positions checkbox depending on if there is peak data
        self.setpeaks_window.inital_peak_info[0] = np.sort(self.setpeaks_window.inital_peak_info[0])
        print(self.setpeaks_window.inital_peak_info)
        if np.any(self.setpeaks_window.inital_peak_info[0]) is None:
            self.show_peakposotions_button.setDisabled(True)
        else:
            #sets the inital peak info to the values calculated in the find peaks window
            self.inital_peak_info = self.setpeaks_window.inital_peak_info
            self.show_peakposotions_button.setDisabled(False)
            #draws the data to the graph
            self.graph_data()
            #corrects for spanning in data indicies
            for i in range(self.peak_no):
                self.inital_peak_info[0][i] = self.inital_peak_info[0][i] + self.sind[0]
            #adds data to the table
            self.itemise_table_from_find_peaks()
        #fit peaks is only enable if there is a background or a guess
        self.enable_fpeaks()
        return

    def open_data_range_window(self):
        """function that allows the user to select a range of data to be graphed"""
        #calls data range function, which adds the spanning window and returns the data range
        self.data_range_window = XPS_DataRangeWindow()
        self.data_range_window.data = self.data
        self.data_range_window.sind = [0,len(self.data[:,0])]
        self.data_range_window.graph_data_basic()
        self.data_range_window.exec_()
        self.sind = self.data_range_window.sind
        self.add_table_data_range()
        return       

    def add_table_columns_from_peak_no(self):
        """clears columns already in the table, and sets to the correct numper of columns with the peak no as their header"""
        self.peak_table.setHorizontalHeaderItem(0,QtWidgets.QTableWidgetItem("Peak {no}".format(no = 0+1)))
        for i in range(1,self.peak_table.columnCount()):
            self.peak_table.removeColumn(1)
        for i in range(1,self.peak_no):
            self.peak_table.insertColumn(i)
            self.peak_table.setHorizontalHeaderItem(i,QtWidgets.QTableWidgetItem("Peak {no}".format(no = i+1)))
        return

    def add_table_peak_positions(self):
        """adds the centre position of each peak to the table"""
        for i in range(self.peak_no):
            item_centre = self.cell(str(self.data[self.inital_peak_info[0][i],0]))
            self.peak_table.setItem(2,i,QtWidgets.QTableWidgetItem(item_centre))
    
    def add_table_data_range(self):
        """adds the data range to the table"""
        for i in range(self.peak_no):
            item_dataL = self.cell(str(self.data[self.sind[0],0]))
            item_dataR = self.cell(str(self.data[self.sind[1]-1,0]))
            self.peak_table.setItem(0,i,QtWidgets.QTableWidgetItem(item_dataL))
            self.peak_table.setItem(1,i,QtWidgets.QTableWidgetItem(item_dataR))
        return
    
    def add_table_fitting_data(self):
        """adds peak data to table """
        for i in range(self.peak_no):
            item_height = self.cell(str(np.round(self.fit_results.params["v{}height".format(i)].value,2)))
            self.peak_table.setItem(3,i,QtWidgets.QTableWidgetItem(item_height))
            item_sigma = self.cell(str(np.round(self.fit_results.params["v{}sigma".format(i)].value,6)))
            self.peak_table.setItem(4,i,QtWidgets.QTableWidgetItem(item_sigma))
            item_gamma = self.cell(str(np.round(self.fit_results.params["v{}gamma".format(i)].value,6)))
            self.peak_table.setItem(5,i,QtWidgets.QTableWidgetItem(item_gamma))
            item_fwhm = self.cell(str(np.round(self.fit_results.params["v{}fwhm".format(i)].value,6)))
            self.peak_table.setItem(6,i,QtWidgets.QTableWidgetItem(item_fwhm))

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
        self.graph_data()
        return
   
    def show_background(self):
        """graphs data if checkbox is clicked"""
        self.graph_data()
        return

    def show_fitting(self):
        """graphs data if checkbox is clicked"""
        self.graph_data()
        return

    def show_fitted_peaks(self):
        """graphs data if checkbox is clicked"""
        self.graph_data()
        return
    
    def initiate_fitting(self):
        """runs the fitting algorithm and set a number of buttons to be enabled or disabled depending on the outcome"""
        for i in range(self.peak_no):
            self.inital_peak_info[0][i] = self.inital_peak_info[0][i] - self.sind[0]
        self.gind[0],self.gind[1] = self.sind[0],self.sind[1]
        #runs the peak fitting algorithm
        self.fit_results = peak_fitting.fit_peaks(self,self.data[int(self.gind[0]):int(self.gind[1]),:],self.peak_no,self.inital_peak_info,self.background[int(self.gind[0]):int(self.gind[1])])
        for i in range(self.peak_no):
            self.inital_peak_info[0][i] = self.inital_peak_info[0][i] + self.sind[0]
        #if fitting has been completed successfully the will enable the checkbox to show the fitting curve and fitted peaks, and the save peak data button
        if self.fit_results is not None:
            self.show_fittingcurve_button.setEnabled(True)
            self.show_fittedpeaks_button.setEnabled(True)
            self.save_button.setEnabled(True)
            self.display_fitting()
            self.graph_data()

        return

    def display_fitting(self):
        """prints out fit report """
        self.display_fit_report()
        self.add_table_fitting_data()
        return

    def enable_fpeaks(self):
        """checks if there is a background or an inital guess and enables the fit peaks button if there is"""
        if self.background is not None and self.inital_peak_info is not None:
            self.fit_button.setDisabled(False)
        else:
            self.fit_button.setDisabled(True)
    
    
    def save_report(self):
        """saves infromation from the peak model into a JSON file, as provided by lmfit """
        file_name, _ = QtWidgets.QFileDialog.getSaveFileName(self,"Save JSON File","", "JSON file (*.json)")
        lm.model.save_model(self.model,file_name)
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
