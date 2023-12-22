from matplotlib.pyplot import text
import lmfit as lm
import numpy as np
from PyQt5 import uic,QtWidgets, QtGui
from PyQt5.QtCore import Qt
from pathlib import Path
import ast
import os
import json
tiny = 1.0e-15

class XPS_ConstraintsWindow(QtWidgets.QDialog):
    def __init__(self):
        super(XPS_ConstraintsWindow,self).__init__()
        self.inital_peak_positions = None
        self.pars = None
        self.prefix = None
        self.peak_no = 1
        self.constraints = None
        self.initUI()
        
    def initUI(self):
        path = Path("UIWidgets")
        filename = path / "XPS_constraints.ui"
        uic.loadUi(filename,self)
        self.constraint_parameters_combobox = self.findChild(QtWidgets.QComboBox,"comboBox_Parameters")        
        self.constraint_parameters_combobox.currentTextChanged.connect(self.update_table)

        self.constraint_peaks_combobox = self.findChild(QtWidgets.QComboBox,"comboBox_Peaks")
        
        self.constraint_type_combobox = self.findChild(QtWidgets.QComboBox,"comboBox_TypeSet")
        self.constraint_type_combobox.currentTextChanged.connect(self.set_value_combo_options)
        
        self.constraint_value_combobox = self.findChild(QtWidgets.QComboBox,"comboBox_Values")
        
        self.set_lineedit = self.findChild(QtWidgets.QLineEdit,"lineEdit_SetValue")
        
        self.clear_button = self.findChild(QtWidgets.QPushButton,"button_ClearParameters")
        self.clear_button.clicked.connect(self.clear_parameters)
        
        self.add_con_button = self.findChild(QtWidgets.QPushButton,"button_AddConstraint")
        self.add_con_button.clicked.connect(self.add_constraint)
        
        self.save_button = self.findChild(QtWidgets.QPushButton,"button_SaveCons")
        self.save_button.clicked.connect(self.save_constraints)
        
        self.save_button = self.findChild(QtWidgets.QPushButton,"button_LoadCons")
        self.save_button.clicked.connect(self.load_constraints)
        
        self.done_button = self.findChild(QtWidgets.QPushButton,"button_Done")
        self.done_button.clicked.connect(self.close)
        
        self.cons_table = self.findChild(QtWidgets.QTableWidget,"tableWidget_Cons")
        self.set_value_combo_options()
        return
   
 
    def update_display_combobox(self):
        self.constraint_parameters_combobox.clear()
        if self.pars is not None:
            keys = self.pars.keys()
        else:
            return
        #remove prefix from keys
        keys = [key[len(self.prefix[0]):] for key in keys]
        items_list = list(set(keys))
        self.constraint_parameters_combobox.addItems(items_list)
        return
    
    def update_table(self):
        param = self.constraint_parameters_combobox.currentText()
        if param == "":
            return
        
        key = list(self.constraints.keys())
        keys = list(self.constraints[key[0]].keys())

        self.cons_table.setRowCount(0)       
        for i in range(self.peak_no):
            self.cons_table.insertRow(i)
            self.cons_table.setVerticalHeaderItem(i,QtWidgets.QTableWidgetItem("Peak {} ".format(i+1) + self.prefix[i]))
            keys = list(self.constraints[self.prefix[i]+param.lower()].keys())
            for j in range(len(keys)):
                if keys[j] == "Vary":
                    item = self.cell(str(self.constraints[self.prefix[i]+param.lower()][keys[j]]))
                if keys[j] == "value" and param == "Center":
                    item = self.center_value_cell(str(self.constraints[self.prefix[i]+param.lower()][keys[j]]))
                else:
                    item = self.cell(str(self.constraints[self.prefix[i]+param.lower()][keys[j]]))
                self.cons_table.setItem(i,j,item)
        self.allow_add_button()
        return
    
    def cell(self,var):
        item = QtWidgets.QTableWidgetItem()
        item.setText(var)
        item.setFlags(item.flags() ^ Qt.ItemIsEditable)
        return item   
     
    def center_value_cell(self,var):
        item = QtWidgets.QTableWidgetItem()
        item.setText(var)
        item.setFlags(item.flags() ^ Qt.ItemIsEditable)
        item.setBackground(Qt.grey)
        return item
        
    def set_peak_number_options(self):
        """Function to set the peak numbers in the combo boxes to the correct number of peaks """
        self.constraint_peaks_combobox.clear()
        if self.pars == None:
            return
        self.constraint_peaks_combobox.addItem("All Peaks")
        for i in range(1,self.peak_no+1):
            self.constraint_peaks_combobox.addItem("Peak {} ".format(i) + "{}".format(self.prefix[i-1]))
        return
       
    
    def set_value_combo_options(self):
        """Function that sets the options for the value combo box"""
        self.constraint_value_combobox.clear()
        if self.constraint_type_combobox.currentText() == "Min" or self.constraint_type_combobox.currentText() == "Max":
            self.constraint_value_combobox.addItems(["% from Value","Exactly"])
            self.set_lineedit.setEnabled(True)
        elif self.constraint_type_combobox.currentText() == "Value":
            self.constraint_value_combobox.addItems(["Exactly"])
            self.set_lineedit.setEnabled(True)
        elif self.constraint_type_combobox.currentText() == "Expr":
            self.set_lineedit.setEnabled(True)
            self.constraint_value_combobox.addItems(["As Expressed"])
        elif self.constraint_type_combobox.currentText() == "Vary":
            self.set_lineedit.setEnabled(False)
            self.constraint_value_combobox.addItems(["True","False"])
        self.allow_add_button()
        return
    
    def add_constraint(self):
        """Function that adds a constraint to the constraints dictionary"""
        if self.set_lineedit.text() == "" and  self.constraint_type_combobox.currentText() != "Vary":
            print("No value entered")
            return
        
        if self.constraint_type_combobox.currentText() == "Expr":
            if self.is_valid_python(self.set_lineedit.text()) == False:
                print("Invalid Python Expression")
                return
            else:
                print("Valid Python Expression")
                try:
                    self.add_expression()
                except:
                    print("unable to add: ",self.set_lineedit.text())
                    print(self.constraint_peaks_combobox.currentText())
                    print(self.constraint_peaks_combobox.currentText()[-len(self.prefix[0]):])
                    return
                    
            
        elif self.constraint_type_combobox.currentText() == "Min" or self.constraint_type_combobox.currentText() == "Max" or self.constraint_type_combobox.currentText() == "Value":
            if self.check_number(self.set_lineedit.text()) == False:
                print("Not a number")
                return
            else:
                self.add_number()
        else:
            self.add_vary()
        self.update_table()
        return
    
    def check_number(self,text):
        """Function that checks if the input is a number"""
        if text[0] == "-":
            text = text[1:]
        else:
            pass
        if text.replace('.','',1).isdigit() == True:
            return True
        else:
            return False
        
    def is_valid_python(self,code):
        try:
            ast.parse(code)
        except SyntaxError:           
            return False
        return True
    
    def add_expression(self):
        """Function that adds an expression to the constraints dictionary"""
        if self.constraint_peaks_combobox.currentText() == "All Peaks":
            for i in range(self.peak_no):
                if str(self.set_lineedit.text()) == "None":
                    self.constraints[self.prefix[i]+self.constraint_parameters_combobox.currentText().lower()]["Expr"] = None
                else:                    
                    self.constraints[self.prefix[i]+self.constraint_parameters_combobox.currentText().lower()]["Expr"] = str(self.set_lineedit.text())
        else:
            if str(self.set_lineedit.text()) == "None":
                self.constraints[self.constraint_peaks_combobox.currentText()[-len(self.prefix[0]):]+self.constraint_parameters_combobox.currentText().lower()]["Expr"] = None
            else:
                self.constraints[self.constraint_peaks_combobox.currentText()[-len(self.prefix[0]):]+self.constraint_parameters_combobox.currentText().lower()]["Expr"] = str(self.set_lineedit.text())
        return
    
    def add_number(self):
        """Function that adds a number to the constraints dictionary"""
        if self.constraint_value_combobox.currentText() == "Exactly":
            if self.constraint_peaks_combobox.currentText() == "All Peaks":
                for i in range(self.peak_no):
                    self.constraints[self.prefix[i]+self.constraint_parameters_combobox.currentText().lower()][self.constraint_type_combobox.currentText()] = float(self.set_lineedit.text())
            else:
                self.constraints[self.constraint_peaks_combobox.currentText()[-len(self.prefix[0]):]+self.constraint_parameters_combobox.currentText().lower()][self.constraint_type_combobox.currentText()] = float(self.set_lineedit.text())
        else:
            if self.constraint_peaks_combobox.currentText() == "All Peaks":
                for i in range(self.peak_no):
                    value = self.constraints[self.prefix[i]+self.constraint_parameters_combobox.currentText().lower()]["Value"]
                    if self.constraint_type_combobox.currentText() == "Min":
                        self.constraints[self.prefix[i]+self.constraint_parameters_combobox.currentText().lower()][self.constraint_type_combobox.currentText()] = float(value)*(1-float(self.set_lineedit.text())/100)
                    else:
                        self.constraints[self.prefix[i]+self.constraint_parameters_combobox.currentText().lower()][self.constraint_type_combobox.currentText()] = float(value)*(1+float(self.set_lineedit.text())/100)  
            else:
                value = self.constraints[self.constraint_peaks_combobox.currentText()[-len(self.prefix[0]):]+self.constraint_parameters_combobox.currentText().lower()]["Value"]
                if self.constraint_type_combobox.currentText() == "Min":                   
                    self.constraints[self.constraint_peaks_combobox.currentText()[-len(self.prefix[0]):]+self.constraint_parameters_combobox.currentText().lower()][self.constraint_type_combobox.currentText()] = float(value)*(1-float(self.set_lineedit.text())/100)
                else:
                    self.constraints[self.constraint_peaks_combobox.currentText()[-len(self.prefix[0]):]+self.constraint_parameters_combobox.currentText().lower()][self.constraint_type_combobox.currentText()] = float(value)*(1+float(self.set_lineedit.text())/100)       
        return
    
    def select_setting(self):
        if self.constraint_type_combobox.currentText() == "Expr":
            self.constraint_value_combobox.setEnabled(False)
            self.set_lineedit.setEnabled(True)
    
    def add_vary(self):
        """Function that adds a vary to the constraints dictionary"""
        if self.constraint_peaks_combobox.currentText() == "All Peaks":
            for i in range(self.peak_no):
                self.constraints[self.prefix[i]+self.constraint_parameters_combobox.currentText().lower()]["Vary"] = bool(self.constraint_value_combobox.currentText())
        else:
            self.constraints[self.constraint_peaks_combobox.currentText()[len(self.prefix[0]):]+self.constraint_parameters_combobox.currentText().lower()]["Vary"] = bool(self.constraint_value_combobox.currentText())
        return
    
    def clear_parameters(self):
        self.allow_add_button()
        self.update_table()
        return
     
    def allow_add_button(self):
        """Function that enables the add button"""
        if self.constraint_parameters_combobox.currentText() == "center" and self.constraint_type_combobox.currentText() == "Value":
            self.add_con_button.setEnabled(False)
        else:
            self.add_con_button.setEnabled(True)
        return
    
    def save_constraints(self):
        """Function that saves the constraints dictionary to a json file"""
        file_name = QtWidgets.QFileDialog.getSaveFileName(self, 'Save File', os.getcwd(), "JSON (*.json)")
        if file_name[0] != "":
            with open(file_name[0], 'w') as fp:
                json.dump(self.constraints, fp)
        return
    
    def load_constraints(self):
        """Function that loads a constraints dictionary from a json file"""
        file_name = QtWidgets.QFileDialog.getOpenFileName(self, 'Open File', os.getcwd(), "JSON (*.json)")
        if file_name[0] != "":
            with open(file_name[0], 'r') as fp:
                self.constraints = json.load(fp)
            self.update_table()
        return
    
    def done1(self):
        self.close()
        return
    
class XPS_SetPeakModel(QtWidgets.QDialog):
    def __init__(self,peak_no,prefix):
        super(XPS_SetPeakModel,self).__init__()
        self.peak_no = peak_no
        self.prefix = prefix
        self.set_model = []
        self.initUI()
    
    def initUI(self):
        self.setWindowTitle("Set Peak Model")
        path = Path("UIWidgets")
        filename = path / "XPS_model.ui"
        uic.loadUi(filename,self)
        self.model_combobox = self.findChild(QtWidgets.QComboBox,"comboBox_Model")
        self.peak_combobox = self.findChild(QtWidgets.QComboBox,"comboBox_Peaks")        
        self.model_table = self.findChild(QtWidgets.QTableWidget,"tableWidget_Model")
        
        self.add_button = self.findChild(QtWidgets.QPushButton,"button_AddModel")
        self.add_button.clicked.connect(self.add_model)
        self.done_button = self.findChild(QtWidgets.QPushButton,"button_Done")
        self.done_button.clicked.connect(self.done1)        
        self.initiate_table_and_peak_box()
    
    def initiate_table_and_peak_box(self):
        """Function that initializes the model table"""
        self.model_table.setRowCount(self.peak_no)
        peak_headers = ["peak_{}".format(i) for i in range(1,self.peak_no+1)]
        self.peak_combobox.clear()
        self.peak_combobox.addItems(peak_headers)
        self.model_table.setColumnCount(1) 
        self.model_table.setVerticalHeaderLabels(peak_headers)
        self.model_table.setHorizontalHeaderLabels(["Model"])
        self.model_table.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.Stretch)
        self.model_table.verticalHeader().setSectionResizeMode(QtWidgets.QHeaderView.Stretch)
        self.model_table.setEditTriggers(QtWidgets.QTableWidget.NoEditTriggers)
        return
    
    def add_model(self):
        """Function that adds a model to the model table"""
        row = self.peak_combobox.currentIndex()
        self.model_table.setItem(row,0,QtWidgets.QTableWidgetItem(self.model_combobox.currentText()))
        self.set_model = self.set_model_list()
        return

    def set_model_list(self):
        set_model = []
        for i in range(self.peak_no):
            if self.model_table.item(i,0) is None:
                pass
            else:
                set_model.append(self.model_table.item(i,0).text())
        return set_model


    def done1(self):
        self.close()
        return
    
class custom_models:
        
    def reversed_DoniachSunjic(x, amplitude=1.0, center=0, sigma=1.0, gamma=0.0):
        """Return a Doniach Sunjic asymmetric lineshape.

        doniach(x, amplitude, center, sigma, gamma) =
            amplitude / sigma^(1-gamma) *
            cos(pi*gamma/2 + (1-gamma) arctan((x-center)/sigma) /
            (sigma**2 + (x-center)**2)**[(1-gamma)/2]

        For example used in photo-emission; see
        http://www.casaxps.com/help_manual/line_shapes.htm for more information.

        """
        arg = -(x-center)/max(tiny, sigma)
        gm1 = (1.0 - gamma)
        scale = amplitude/max(tiny, (sigma**gm1))
        return scale*np.cos(np.pi*gamma/2 + gm1*np.arctan(arg))/(1 + arg**2)**(gm1/2)

    
    def DoniachSunjic_Guassian_Conv(x, damplitude=1.0, dcenter=0, dsigma=1.0, dgamma=0.0,gamplitude=1.0, gcenter=0.0, gsigma=1.0):
        arg = (x-dcenter)/max(tiny, dsigma)
        gm1 = (1.0 - dgamma)
        scale = damplitude/max(tiny, (dsigma**gm1))
        
        D = scale*np.cos(np.pi*dgamma/2 + gm1*np.arctan(arg))/(1 + arg**2)**(gm1/2)
        G = ((gamplitude/(max(tiny, np.pi*gsigma)))*np.exp(-(1.0*x-gcenter)**2 / max(tiny, (2*gsigma**2))))  
        C = np.convolve(D,G,'same')
        return C
    
    def reversed_DoniachSunjic_Guassian_Conv(x, damplitude=1.0, dcenter=0, dsigma=1.0, dgamma=0.0,gamplitude=1.0, gcenter=0.0, gsigma=1.0):
        arg = -(x-dcenter)/max(tiny, dsigma)
        gm1 = (1.0 - dgamma)
        scale = damplitude/max(tiny, (dsigma**gm1))
        
        D = scale*np.cos(np.pi*dgamma/2 + gm1*np.arctan(arg))/(1 + arg**2)**(gm1/2)
        G = ((gamplitude/(max(tiny, np.pi*gsigma)))*np.exp(-(1.0*x-gcenter)**2 / max(tiny, (2*gsigma**2))))  
        C = np.convolve(D,G,'same')
        return C
    
class peak_fitting(custom_models):
    """Peak fitting class provides functions for finding the best model to fit to the XPS data"""
    
    def set_individual_peaks(self,peak_no,set_model_list):
        peaks = np.empty(peak_no,dtype = "object")
        ext = 0
        prefix_list = []
        for i in range(peak_no):
            peaks[i],prefix = self.get_model(i+ext,set_model_list[i],peak_no)
            ext += 1
            prefix_list.extend(prefix)
        return peaks,prefix_list

    def get_model(self,i,model_name,peak_no):
        if model_name == "Voigt Peak":
            prefix = "V" + str(i)
            model = lm.models.VoigtModel(prefix = prefix)
        elif model_name == "Lorentzian Peak":
            prefix = "L" + str(i)
            model = lm.models.LorentzianModel(prefix = prefix)  
        elif model_name == "Doniach-Sunjic Peak":
            prefix = "D" + str(i)
            model = lm.models.DoniachModel(prefix = prefix)
        elif model_name == "reversed Doniach-Sunjic Peak":
            prefix = "rD" + str(i)
            print("reversed Doniach-Sunjic Peak")            
            model = lm.models.Model(custom_models.DoniachSunjic_Guassian_Conv,prefix = prefix)
            print(model)
        elif model_name == "Convolved Doniach-Sunjic and Guassian":
            prefix = "C" + str(i)
            model = lm.models.Model(custom_models.DoniachSunjic_Guassian_Conv,prefix = prefix)
        elif model_name == "reversed Convolved Doniach-Sunjic and Guassian":
            prefix = "rC" + str(i)
            model = lm.models.Model(custom_models.reversed_DoniachSunjic_Guassian_Conv,prefix = prefix)
        else:
            print("model: ",model_name, " is not recognised")
        return model,prefix
    
    def VoigtModel_lmfit(self,peak_no,center_pos):
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
    
    def LorentzianModel_lmfit(self,peak_no,center_pos):
        """Function that creates an object array, each representing one peak"""
        #initializes an empty array, before populating each element with a VoightModel object
        #Each peak is labeled in the fit report by the prefix v{i}
        peaks = np.empty(peak_no,dtype = "object")
        model_prefix = ["L{}".format(i) for i in range(peak_no)]   
        for i in range(peak_no):
            peaks[i] = lm.models.LorentzianModel(prefix = model_prefix[i])
        model = np.sum(peaks)
        pars = model.make_params()
        for i in range(len(peaks)):
            pars[model_prefix[i] + "amplitude"].set(value = 10,min = 1)
            pars[model_prefix[i]+"center"].set(center_pos[i],min = center_pos[i]-0.1,max = center_pos[i]+0.1)
            pars[model_prefix[i]+"sigma"].set(value = 0.1,min = 0.01,max = 0.5)
        return model,pars,model_prefix

    def DoniachModel_lmfit(self,peak_no,center_pos):
        """Function that creates an object array, each representing one peak"""
        #initializes an empty array, before populating each element with a VoightModel object
        #Each peak is labeled in the fit report by the prefix v{i}
        peaks = np.empty(peak_no,dtype = "object")
        model_prefix = ["D{}".format(i) for i in range(peak_no)]
        for i in range(peak_no):
            peaks[i] = lm.models.DoniachModel(prefix = model_prefix[i])
        model = np.sum(peaks)
        pars = model.make_params()
        for i in range(len(peaks)):
            pars[model_prefix[i] + "amplitude"].set(value = 10,min = 1)
            pars[model_prefix[i] + "center"].set(center_pos[i],min = center_pos[i]-0.1,max = center_pos[i]+0.1)
            pars[model_prefix[i] + "sigma"].set(value = 0.1,min = 0.01,max = 0.5)
            pars[model_prefix[i] + "gamma"].set(value = 0.1,min = 0.01,max = 0.5,expr = None,vary = True)
        return model,pars,model_prefix

    def reversed_DoniachModel_lmfit(self,peak_no,center_pos):
        """Function that creates an object array, each representing one peak"""
        #initializes an empty array, before populating each element with a VoightModel object
        #Each peak is labeled in the fit report by the prefix v{i}
        peaks = np.empty(peak_no,dtype = "object")
        model_prefix = ["rD{}".format(i) for i in range(peak_no)]
        for i in range(peak_no):
            peaks[i] = lm.models.Model(custom_models.reversed_DoniachSunjic,independent_vars=['x'],prefix = model_prefix[i])
        model = np.sum(peaks)
        pars = model.make_params()
        for i in range(len(peaks)):
            pars[model_prefix[i] + "amplitude"].set(value = 10,min = 1)
            pars[model_prefix[i] + "center"].set(center_pos[i],min = center_pos[i]-0.1,max = center_pos[i]+0.1)
            pars[model_prefix[i] + "sigma"].set(value = 0.1,min = 0.01,max = 0.5)
            pars[model_prefix[i] + "gamma"].set(value = 0.1,min = 0.01,max = 0.5,expr = None,vary = True)
        return model,pars,model_prefix
    
    def Convolved_DoniachGuassian_lmfit(self,peak_no,center_pos):
        """Function that creates an object array, each representing one peak"""
        #initializes an empty array, before populating each element with a VoightModel object
        #Each peak is labeled in the fit report by the prefix v{i}
        peaks = np.empty(peak_no,dtype = "object")
        model_prefix = ["C{}".format(i) for i in range(peak_no)]
        for i in range(peak_no):
            peaks[i] = lm.models.Model(custom_models.DoniachSunjic_Guassian_Conv,prefix = model_prefix[i])
        model = np.sum(peaks)
        pars = model.make_params()
        for i in range(len(peaks)):
            pars[model_prefix[i] + "damplitude"].set(value = 10,min = 1)
            pars[model_prefix[i] + "dcenter"].set(center_pos[i],min = center_pos[i]-0.1,max = center_pos[i]+0.1)
            pars[model_prefix[i] + "dsigma"].set(value = 0.1,min = 0.01,max = 0.5)
            pars[model_prefix[i] + "dgamma"].set(value = 0.1,min = 0.01,max = 0.5,expr = None,vary = True)
            pars[model_prefix[i] + "gamplitude"].set(value = 10,min = 1,expr = model_prefix[i] + "damplitude",vary = False)
            pars[model_prefix[i] + "gcenter"].set(expr = model_prefix[i] + "dcenter")
            pars[model_prefix[i] + "gsigma"].set(value = 0.1,min = 0.01,max = 0.5)            
        return model,pars,model_prefix

    def reveresed_Convolved_DoniachGuassian_lmfit(self,peak_no,center_pos):
        """Function that creates an object array, each representing one peak"""
        #initializes an empty array, before populating each element with a VoightModel object
        #Each peak is labeled in the fit report by the prefix v{i}
        peaks = np.empty(peak_no,dtype = "object")
        model_prefix = ["rC{}".format(i) for i in range(peak_no)]
        for i in range(peak_no):            
            peaks[i] = lm.models.Model(custom_models.reversed_DoniachSunjic_Guassian_Conv,prefix = model_prefix[i])
        model = np.sum(peaks)
        pars = model.make_params()
        for i in range(len(peaks)):
            pars[model_prefix[i] + "damplitude"].set(value = 10,min = 1)
            pars[model_prefix[i] + "dcenter"].set(center_pos[i],min = center_pos[i]-0.1,max = center_pos[i]+0.1)
            pars[model_prefix[i] + "dsigma"].set(value = 0.1,min = 0.01,max = 0.5)
            pars[model_prefix[i] + "dgamma"].set(value = 0.1,min = 0.01,max = 0.5,expr = None,vary = True)
            pars[model_prefix[i] + "gamplitude"].set(value = 10,min = 1,expr = model_prefix[i] + "damplitude",vary = False)
            pars[model_prefix[i] + "gcenter"].set(expr = model_prefix[i] + "dcenter")
            pars[model_prefix[i] + "gsigma"].set(value = 0.1,min = 0.01,max = 0.5)            
        return model,pars,model_prefix
    
    def convolve(self,dat, kernel):
        """simple convolution of two arrays"""
        npts = min(len(dat), len(kernel))
        pad = np.ones(npts)
        tmp = np.concatenate((pad*dat[0], dat, pad*dat[-1]))
        out = np.convolve(tmp, kernel, mode='valid')
        noff = int((len(out) - npts) / 2)
        return (out[noff:])[:npts]
    
    def set_constraints_variable(self,pars):        
        """Function to set the form of constraints dictionary"""
        if pars == None:
            return None
        constraints = {}
        for param in pars:
            constraints[param] = {"Value":pars[param].value,"Min":pars[param].min,"Max":pars[param].max,"Expr":pars[param].expr,"Vary":pars[param].vary}
        return constraints
    
    def convert_constraints_to_new_model(self):
        """Function that converts the usable parameters from the old model to the new model"""
        #Creates a dictionary of the parameters from the old model
        old_cons = self.constraints
        #Creates a dictionary of the parameters from the new model
        new_cons = self.set_constraints_variable(self.pars)
        #If the parameter is in the new model it is added to the usable parameters dictionary
        for param in old_cons.keys():
            if any(param[2:] in key for key in new_cons.keys()):
                for i in range(self.peak_no):                    
                    new_cons[self.model_prefix[i] + param[2:]]["Value"] = old_cons[param[0] + str(i) + param[2:]]["Value"]
                    new_cons[self.model_prefix[i] + param[2:]]["Min"] = old_cons[param[0] + str(i) + param[2:]]["Min"]
                    new_cons[self.model_prefix[i] + param[2:]]["Max"] = old_cons[param[0] + str(i) + param[2:]]["Max"]
        self.constraints = new_cons
        return
    
    def keep_parameters_fixed(self,pars,contraints):
        """Function to keep the parameters fixed for the peak fitting"""
        #Creates the model and parameters for the peak fitting from scratch if there is none
        if pars is None:
            pass
        elif contraints is None:
            pass
        elif contraints == {}:
            pass
        elif all(contraints.keys()) != all(pars.keys()):
            pass
        else:
            for param in contraints:
                pars[param].set(value = contraints[param]["Value"],min = contraints[param]["Min"],max = contraints[param]["Max"],expr = contraints[param]["Expr"],vary = contraints[param]["Vary"])
        return pars
    
    def generate_model_with_guess(self,x,y,peaks,peak_no):
        """Function to generate the model and parameters for the peak fitting"""
        #Creates the model and parameters for the peak fitting
        #sets the variable model to contain the combine peak information
        model = peaks[0]
        for i in range(1,peak_no):
            model = model + peaks[i]
       #sets the variable pars to contain the parameters for the peak fitting
        pars = model.make_params()
        pars = model.guess(y,x = x)
        return model,pars
    
    def inital_params_voigt(self,pars,center_pos):
        """Function to generate the intial starting parameters for a voigt model"""
        #Creates the model and parameters for the peak fitting
        #sets the variable pars to contain the parameters for the peak fitting
        for i in range(self.peak_no):
            pars["V{}center".format(i)].set(value = center_pos[i])
            pars["V{}sigma".format(i)].set(value = 0.1, min = 0, max = 1)
            pars["V{}amplitude".format(i)].set(value = 1.0,min = 0)
            pars["V{}gamma".format(i)].set(value = 0.1, min = 0, max = 1,expr = "None",vary = True)
        return

    def inital_params_lorenz(self,pars,center_pos,peak_no):
        """Function to generate the intial starting parameters for a voigt model"""
        #Creates the model and parameters for the peak fitting
        #sets the variable pars to contain the parameters for the peak fitting
        for i in range(peak_no):
            pars["L{}center".format(i)].set(value = center_pos[i])
            pars["L{}sigma".format(i)].set(value = 0.1, min = 0, max = 1)
            pars["L{}amplitude".format(i)].set(value = 1.0,min = 0)
        return
    
    def generate_model_with_params(self,x,peaks,peak_no,inital_peak_positions):
        """Function to generate the model and parameters for the peak fitting"""
        #Creates the model and parameters for the peak fitting
        pars = peaks[0].make_params(center = x[inital_peak_positions[0]])      
        #sets the variable model to contain the combine peak information
        model = peaks[0]
        for i in range(1,peak_no):
            pars.update(peaks[i].make_params(center = x[inital_peak_positions[i]]))
            model = model + peaks[i]
        return model,pars

    def generate_model(self,peaks,peak_no):
        """Function to generate the model and parameters for the peak fitting"""
        #Creates the model and parameters for the peak fitting   
        #sets the variable model to contain the combine peak information
        model = peaks[0]
        for i in range(1,peak_no):
            model = model + peaks[i]
        return model     
 
    def set_constraint_parameters(self):
        """Functions that sets the parameters for the constraints tab"""
        for key_param in self.constraints:
            self.pars[key_param].set(value = self.constraints[key_param]["Value"], min = self.constraints[key_param]["Min"], max = self.constraints[key_param]["Max"],expr = self.constraints[key_param]["Expr"],vary = self.constraints[key_param]["Vary"])        
        return
    
    def fit_peaks(self,x,y,pars,model,method):
        #fits the peak model to the XPS data
        print(method)
        fit_results = model.fit(y,pars,x=x,method=method,nan_policy='omit')  
        #print(fit_results.fit_report())
        return fit_results
  
    
    def fit_peaks_sequential(self,x,y,pars,model):
        """Function to fit the peak model to the XPS data, by squentially dueucing each paramter"""
        for i in range(10):
            for param in pars:
                pars[param].vary = False
            for param in pars:
                pars[param].vary = True
                if param[2:] == "height" or param[2:] == "fwhm":
                    pass
                else:
                    print("now doing : ",param,"inital value: ",pars[param].value, end = "final value: ")                
                    fit_results = model.fit(y,pars,x=x,nan_policy='omit')
                    pars[param].value = fit_results.best_values[param]
                    print(pars[param].value)
        return fit_results
    
   
    
    def fit_peaks_findgamma(self,x,y,pars,model):
        gc = 0.5
        s = 2
        step = gc/s
        for i in range(15):
            #fits the peak model to the XPS data
            for j in range(self.peak_no):
                pars[self.model_prefix[j]+"gamma"].expr = str(gc) + "*" + self.model_prefix[j] + "sigma"
            fit_results_centre = model.fit(y,pars,x=x,nan_policy='omit')
            s=2
            g1,g2 = gc - step, gc + step
            print(i,step,g1,g2,gc)            
            for j in range(self.peak_no):
                pars[self.model_prefix[j]+"gamma"].expr = str(g1) + "*" + self.model_prefix[j] + "sigma"
            fit_results_left = model.fit(y,pars,x=x,nan_policy='omit')
            for j in range(self.peak_no):
                pars[self.model_prefix[j]+"gamma"].expr = str(g2) + "*" + self.model_prefix[j] + "sigma"
            fit_results_right = model.fit(y,pars,x=x,nan_policy='omit')
            if fit_results_centre.redchi < fit_results_left.redchi and fit_results_centre.redchi < fit_results_right.redchi:
                gc = gc
                s+=1
            elif fit_results_left.redchi < fit_results_centre.redchi and fit_results_left.redchi < fit_results_right.redchi:
                gc = g1
            else:
                gc = g2
            step = step/s
        fit_results_centre = model.fit(y,pars,x=x,nan_policy='omit')    
        print("\nFitting completed!\n.")
        return fit_results_centre
    
    def plot_fitting(self,x,fitted_array):
        """Function to display the fitting data, including the fit report following fitting calculation"""
        #plots the data and the fitted curve, along with printing the fit report
        
        self.MplWidget.canvas.axes.plot(x,fitted_array)
        self.MplWidget.canvas.draw()
        return

    def display_fit_report(self):
        print(self.fit_results.fit_report())
        return

    def plot_peak_model(self,x,peaks,peak_no):
        for i in range(peak_no):
            self.MplWidget.canvas.axes.fill(x,peaks[i],alpha = 0.5)
        self.MplWidget.canvas.draw()
        return

    def plot_peak_model_cropped(self,x,peaks,peak_no,sind):
        for i in range(peak_no):
            self.MplWidget.canvas.axes.fill(x[sind[0]:sind[1]],peaks[i][sind[0]:sind[1]],alpha = 0.5)
        self.MplWidget.canvas.draw()        
        return



