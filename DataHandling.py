import numpy as np
import matplotlib.pyplot as plt
import vamas

import scipy.signal as sps
from pathlib import Path
from PyQt5 import uic,QtWidgets
from PyQt5.QtCore import Qt
from matplotlib.widgets import SpanSelector

class XPS_DataRangeWindow(QtWidgets.QDialog):
    def __init__(self):
        super(XPS_DataRangeWindow,self).__init__()
        self.data = None
        self.sind = np.zeros(2)
        self.initUI()
        
    def initUI(self):
        path = Path("UIWidgets")
        filename = path / "XPS_datarange.ui"
        uic.loadUi(filename,self)
        
        #defines the table where the data range is displayed
        self.range_table = self.findChild(QtWidgets.QTableWidget,"tableWidget_Range")
        
        self.range_button = self.findChild(QtWidgets.QPushButton,"button_DataRange")
        self.range_button.clicked.connect(self.select_data_range)
        self.done_button = self.findChild(QtWidgets.QPushButton,"button_Done")
        self.done_button.clicked.connect(self.done1) 
        
    def graph_data_basic(self):
        """Function to plot chosen XPS plots data set"""
        #initalizes variables X & Y for plotting data
        X,Y = self.data[:,0],self.data[:,1]
        self.MplWidget.canvas.axes.clear()
        self.MplWidget.canvas.axes.plot(X,Y)
        #draws
        self.MplWidget.canvas.draw()
        return
    
    def select_data_range(self):
        """function that allows the user to select a range of data to be graphed"""
        #calls data range function, which adds the spanning window and returns the data range
        self.data_range()
        return           
    
    def add_table_data_range(self):
        """adds the data range to the table"""
        item_dataL = self.cell(str(self.data[self.sind[0],0]))
        item_dataR = self.cell(str(self.data[self.sind[1],0]))
        self.range_table.setItem(0,0,QtWidgets.QTableWidgetItem(item_dataL))
        self.range_table.setItem(1,0,QtWidgets.QTableWidgetItem(item_dataR))
        return
    
    def data_range(self):
         self.span = SpanSelector(self.MplWidget.canvas.axes, self.onselect, 'horizontal',
                                        useblit=True, interactive=True, drag_from_anywhere=True, props=dict(alpha=0.5, facecolor='red'))
         return
         
    
    def onselect(self, xmin, xmax):
        idx1,idx2 = (np.abs(self.data[:,0] - xmin)).argmin(),(np.abs(self.data[:,0] - xmax)).argmin()
        self.sind = np.sort(np.array([idx1, idx2]))
        self.add_table_data_range()
        self.range_button.setDisabled(True)
        return    
    
    def cell(self,var):
        item = QtWidgets.QTableWidgetItem()
        item.setText(var)
        item.setFlags(item.flags() ^ Qt.ItemIsEditable)
        return item
        
    def done1(self):
        self.close()
        return     
    
class XPS_PeakFindingWindow(QtWidgets.QDialog):
    def __init__(self):
        super(XPS_PeakFindingWindow,self).__init__()
        self.data = None
        self.selected_parameter = "height"
        self.peak_search_parameters = {"height" : None, "distance" : None, "prominence" : None, "width" : None}
        self.peak_no = 0
        self.peak_text = None
        self.inital_peak_positions = None
        self.initUI()
        
    def initUI(self):
        path = Path("UIWidgets")
        filename = path / "XPS_peaksearch.ui"
        uic.loadUi(filename,self)

        self.parameters_combo = self.findChild(QtWidgets.QComboBox,"comboBox_Parameters")
        self.parameters_combo.currentTextChanged.connect(self.set_selected_parameter)

        self.clear_search = self.findChild(QtWidgets.QPushButton,"button_ClearParameters")
        self.clear_search.clicked.connect(self.clear_peak_search_parameters)

        self.height_label = self.findChild(QtWidgets.QLabel,"Height")
        self.distance_label = self.findChild(QtWidgets.QLabel,"Distance")
        self.prominence_label = self.findChild(QtWidgets.QLabel,"Prominence")
        self.width_label = self.findChild(QtWidgets.QLabel,"Width")

        self.parameters_edit = self.findChild(QtWidgets.QLineEdit,"Edit_Parameters")
        self.parameters_edit.textChanged.connect(self.set_peak_search_parameters)

        self.warning_label = self.findChild(QtWidgets.QLabel,"Warning")
        self.warning_label.setVisible(False)

        self.find_button = self.findChild(QtWidgets.QPushButton,"button_Done")
        self.find_button.clicked.connect(self.done1)

        self.find_button = self.findChild(QtWidgets.QPushButton,"button_FPeaks")
        self.find_button.clicked.connect(self.findpeaks)

        self.peaks_text = self.findChild(QtWidgets.QPlainTextEdit,"plainTextEdit_Peaks")

    def set_selected_parameter(self):
        self.selected_parameter = self.parameters_combo.currentText()
        return

    def set_peak_search_parameters(self):
        if self.parameters_edit.text().isdigit() or self.parameters_edit.text() == None:
            self.peak_search_parameters[self.selected_parameter] = int(self.parameters_edit.text())
            self.warning_label.setVisible(False)
        else:
            self.peak_search_parameters[self.selected_parameter] = self.parameters_edit.text()
            self.warning_label.setVisible(True)
        self.update_serach_labels()
        return

    def update_serach_labels(self):
        self.height_label.setText("Height = " + str(self.peak_search_parameters["height"]))
        self.distance_label.setText("Distance = " + str(self.peak_search_parameters["distance"]))
        self.prominence_label.setText("Prominence = " + str(self.peak_search_parameters["prominence"]))
        self.width_label.setText("Width = " + str(self.peak_search_parameters["width"]))
        return
    
    def clear_peak_search_parameters(self):
        self.parameters_edit.clear()
        self.peak_search_parameters["height"] = None
        self.peak_search_parameters["distance"] = None
        self.peak_search_parameters["prominence"] = None
        self.peak_search_parameters["width"] = None       
        self.update_serach_labels()
        self.peaks_text.setPlainText("")
        self.warning_label.setVisible(False)
        return

    def findpeaks(self):
        if (str(self.peak_search_parameters["height"]).replace('.','',1).isdigit() == False and self.peak_search_parameters["height"] != None) or (str(self.peak_search_parameters["distance"]).replace('.','',1).isdigit() == False and self.peak_search_parameters["distance"] != None) or (str(self.peak_search_parameters["prominence"]).replace('.','',1).isdigit() == False and self.peak_search_parameters["prominence"] != None) or (str(self.peak_search_parameters["width"]).isdigit() == False and self.peak_search_parameters["width"] != None):
            self.warning_label.setVisible(True)
            return
        else:
            self.inital_peak_positions, _ = sps.find_peaks(self.data[:,1], height = float(self.peak_search_parameters["height"]), distance = self.peak_search_parameters["distance"], prominence = self.peak_search_parameters["prominence"], width = self.peak_search_parameters["width"])
            data_handling.set_peak_data(self,self.data,self.peak_search_parameters,self.inital_peak_positions)
            self.peak_no = len(self.inital_peak_positions)

            self.peaks_text.setPlainText(self.peak_text)
        return

    def done1(self):
        self.close()
        return

class XPS_PeakSettingWindow(QtWidgets.QDialog):
    """This class is defines the window for allowing the user to manually set the peak positions"""
    def __init__(self):
        super(XPS_PeakSettingWindow,self).__init__()
        self.data = None
        self.spanned_data = None
        self.peak_no = 1
        self.inital_peak_positions = np.array([]).astype(int)
        self.initUI()
        
    def initUI(self):
        path = Path("UIWidgets")
        filename = path / "XPS_peaksearchman.ui"
        uic.loadUi(filename,self)
        
        self.peaks_combo = self.findChild(QtWidgets.QComboBox,"comboBox_PeakNo")
        self.peaks_combo.currentTextChanged.connect(self.set_peak_no)

        self.done_button = self.findChild(QtWidgets.QPushButton,"button_Done")
        self.done_button.clicked.connect(self.done1)        
        
        self.peak_set_table = self.findChild(QtWidgets.QTableWidget,"tableWidget_SetPeaks")
        self.peak_set_table.itemChanged.connect(self.item_check)
               
    def set_peak_no(self):
        """sets the correct number of rows to be set by the user in the table"""
        current_rows = self.peak_set_table.rowCount()
        self.peak_no = int(self.peaks_combo.currentText())
        self.peak_set_table.setVerticalHeaderItem(0,QtWidgets.QTableWidgetItem("Peak {no}".format(no = 0+1)))
        if self.peak_no > current_rows:
            for i in range(current_rows,self.peak_no):
                self.peak_set_table.insertRow(i)
                self.peak_set_table.setVerticalHeaderItem(i,QtWidgets.QTableWidgetItem("Peak {no}".format(no = i+1)))
        elif self.peak_no < current_rows:
            for i in range(current_rows - self.peak_no):
                self.peak_set_table.removeRow(current_rows - i-1)
        else:
            pass
        self.set_inital_peak_positions_length()
        return

    def item_check(self, item):
        """checks that the user has entered a valid value for the peak position"""
        if item is None:
            pass
        elif item.text().replace(".","",1).isdigit() == True or (item.text()[1:].replace(".","",1).isdigit() == True and item.text()[0] == "-"):
            print(float(item.text()))
            if self.data[0,0] < self.data[-1,0]:
                if (float(item.text()) >= self.data[0,0]) & (float(item.text()) <= self.data[-1,0]):
                    
                    self.set_inital_peak_positions_values()
                else:
                    self.peak_set_table.setItem(item.row(), item.column(), None)
            else:
                if (float(item.text()) <= self.data[0,0]) & (float(item.text()) >= self.data[-1,0]):
                    self.set_inital_peak_positions_values()
                else:
                    self.peak_set_table.setItem(item.row(), item.column(), None)
        else:
            self.peak_set_table.setItem(item.row(), item.column(), None)
        self.graph_data_basic_plus_points()
        return
    
    def set_inital_peak_positions_length(self):
        if len(self.inital_peak_positions) > self.peak_set_table.rowCount():
            self.inital_peak_positions = np.delete(self.inital_peak_positions,np.where(self.inital_peak_positions == 0))
        elif len(self.inital_peak_positions) < self.peak_set_table.rowCount():
            self.inital_peak_positions = np.append(self.inital_peak_positions,np.zeros(self.peak_set_table.rowCount() - len(self.inital_peak_positions)))
        else:
            pass
        return
    
    def set_inital_peak_positions_values(self):
        for i in range(self.peak_set_table.rowCount()):
            if self.peak_set_table.item(i,0) is not None:
                self.inital_peak_positions[i] = np.abs(self.data[:,0] - float(self.peak_set_table.item(i,0).text())).argmin()
            else:
                pass
        return
    
    def graph_data_basic_plus_points(self):
        """Function to plot chosen XPS plots data set"""
        #initalizes variables X & Y for plotting data
        X,Y = self.spanned_data[:,0],self.spanned_data[:,1]
        self.MplWidget.canvas.axes.clear()
        self.MplWidget.canvas.axes.plot(X,Y)
        if np.any(self.inital_peak_positions) is None:
            pass
        else:
            for i in range(len(self.inital_peak_positions)):
                self.MplWidget.canvas.axes.plot(X[int(self.inital_peak_positions[i])],Y[int(self.inital_peak_positions[i])],'x',color = "green")
        #draws
        self.MplWidget.canvas.draw()
        return

    def done1(self):       
        self.close()
        return
    
class data_handling():
    """Data handling class provides functions for manipulating and seperating XPS data for further analysis"""
    def open_vamas_file(filename):
        """Takes the filepath as input and returns an object containing all the vamas file information """
        return vamas.VAMAS(filename)
    
    def select_data_axes(vms, axis):
        """generates data for x & y axes, with selection of which x axis data variant is displayed"""
        #intializes y-axis data
        Y = vms.blocks[0].data[0]
        #initalizes single data array for both x & y values
        data = np.zeros((len(Y),len(Y)))
        #determins which x axis data is being chosen before storing in data, returns 1 if axis is unrecognized
        if axis == "axis":
            X = vms.blocks[0].axis
        elif axis == "kinetic_axis":
            X = vms.blocks[0].kinetic_axis
        elif axis == "binding_axis":
            X = vms.blocks[0].binding_axis
        else:
            print("data format not recognized")
            return 1
        data = np.zeros((len(Y),2))
        data[:,0],data[:,1] = X,Y
        return data

    def graph_data(self):
        """Function to plot chosen XPS plots data set"""
        #initalizes variables X & Y for plotting data
        X,Y = self.data[:,0],self.data[:,1]
        self.MplWidget.canvas.axes.clear()
        self.MplWidget.canvas.axes.plot(X,Y)
        #function includes statments so background is plotted if available
        if self.background is None or self.show_background_button.isChecked() == False:
            pass
        else:
            self.MplWidget.canvas.axes.plot(X,self.background)
        #function includes statments so inital peak positions are plotted if available
        if self.inital_peak_positions is None or self.show_peakposotions_button.isChecked() == False:
            pass
        else:
            for i in range(len(self.inital_peak_positions)):
                self.MplWidget.canvas.axes.plot(X[self.inital_peak_positions[i]],Y[self.inital_peak_positions[i]],'x',color = "green")
        #function includes statments so fitting curve is plotted if available
        if self.fit_results is None or self.show_fittingcurve_button.isChecked() == False:
            pass
        else:
            self.plot_fitting(self.data[:,0],self.fitted_array)
        #function includes statments so fitted peaks are plotted if available
        if self.fit_results is None or self.show_fittedpeaks_button.isChecked() == False:
            pass
        else:
            self.plot_peak_model(self.data[:,0],self.peaks,self.peak_no)#,self.background,self.sind)
        #draws
        self.MplWidget.canvas.draw()
        return

    def graph_data_cropped(self):
        """Function to plot chosen XPS plots data set"""
        #initalizes variables X & Y for plotting data
        X,Y = self.data[:,0],self.data[:,1]
        self.MplWidget.canvas.axes.clear()
        self.MplWidget.canvas.axes.plot(X[self.sind[0]:self.sind[1]],Y[self.sind[0]:self.sind[1]])
        #function includes statments so background is plotted if available
        if self.background is None or self.show_background_button.isChecked() == False:
            pass
        else:
            self.MplWidget.canvas.axes.plot(X[self.sind[0]:self.sind[1]],self.background[self.sind[0]:self.sind[1]])
        #function includes statments so inital peak positions are plotted if available
        if self.inital_peak_positions is None or self.show_peakposotions_button.isChecked() == False:
            pass
        else:
            for i in range(len(self.inital_peak_positions)):
                if self.inital_peak_positions[i] > self.sind[0] and self.inital_peak_positions[i] < self.sind[1]:
                    self.MplWidget.canvas.axes.plot(X[self.inital_peak_positions[i]],Y[self.inital_peak_positions[i]],'x',color = "green")
        #function includes statments so fitting curve is plotted if available
        if self.fit_results is None or self.show_fittingcurve_button.isChecked() == False:
            pass
        else:
            self.plot_fitting(self.data[self.sind[0]:self.sind[1],0],self.fitted_array[self.sind[0]:self.sind[1]])
        #function includes statments so fitted peaks are plotted if available
        if self.fit_results is None or self.show_fittedpeaks_button.isChecked() == False:
            pass
        else:
            self.plot_peak_model_cropped(self.data[:,0],self.peaks,self.peak_no,self.sind)#,self.background,self.sind)
        #draws
        self.MplWidget.canvas.draw()
        return
    
    def set_peak_data(self,data,peaks,peakx):
        """Function that displays the peak data for that detected using the find peaks function """
        self.peak_text = "Number of peaks detected is: ".format(len(peakx))
        self.peak_text = self.peak_text + "\nPeak data is as follows:"
        for i in range(len(peakx)):
            self.peak_text = self.peak_text + "\npeak "+ str(i+1) + ": position = "+ str(data[peakx[i],0]) + ", height = " + str(data[peakx[i],1])
        return
    
    def set_graph_variables(self):
        #stores the fitted data in an array variable for graphing and saving
        self.fitted_array[int(self.sind[0]):int(self.sind[1])] = self.fit_results.best_fit + self.background[int(self.sind[0]):int(self.sind[1])]
        self.peaks = np.empty((self.peak_no,len(self.data[:,0])))
        self.peaks[:] = np.nan
        comps = self.fit_results.eval_components()
        for i in range(self.peak_no):
            if type(self.model_prefix[0]) == str:
                peak = comps[self.model_prefix[i]]+self.background[self.sind[0]:self.sind[1]]
                self.peaks[i][self.sind[0]:self.sind[1]] = peak
            else:
                for peak_name in self.model_prefix[i]:
                    peak = comps[peak_name]+self.background[self.sind[0]:self.sind[1]]
                    self.peaks[i][self.sind[0]:self.sind[1]] = peak
        return
    
    def normalise_data(data):
        """Function to normalise the data to be in a range of 0-1000"""
        #normalises data
        data[:,1] = data[:,1]/max(data[:,1])*1000
        return
    
    def check_vms(vms):
        dirt = np.array(dir(vms))
        for i in range(len(dirt)):
            print(dirt[i])
            print(getattr(vms,dirt[i]))
            print(" ")
            print(" ")
        return
