import numpy as np
import numpy.ma as ma
from PyQt5 import uic,QtWidgets
from pathlib import Path
import sys
import warnings
from matplotlib.widgets import SpanSelector


class xps_BackgroundWindow(QtWidgets.QDialog):
	def __init__(self):
		super(xps_BackgroundWindow,self).__init__()
		self.data = None
		self.sind = np.zeros(2).astype(int)
		self.background_type_1 = "None"
		self.background_range_1 = np.zeros(2).astype(int)
		self.background = None
		self.background_type_2 = "None"
		self.background_range_2 = np.zeros(2).astype(int)
		self.background_2 = None
		self.initUI()
        
	def initUI(self):
		path = Path("UIWidgets")
		filename = path / "XPS_background.ui"
		uic.loadUi(filename,self)
		self.background_combo = self.findChild(QtWidgets.QComboBox,"comboBox_Background")
		self.background_combo_2 = self.findChild(QtWidgets.QComboBox,"comboBox_Background_2")
		
		self.range_button1 = self.findChild(QtWidgets.QPushButton,"button_Range1")
		self.range_button1.clicked.connect(self.select_range1)
		
		self.range_button2 = self.findChild(QtWidgets.QPushButton,"button_Range2")
		self.range_button2.clicked.connect(self.select_range2)
		
		self.background_button = self.findChild(QtWidgets.QPushButton,"button_FBackground")
		self.background_button.clicked.connect(self.get_background)
		
		self.done_button = self.findChild(QtWidgets.QPushButton,"button_Done") 
		self.done_button.clicked.connect(self.done1)
		
		self.show_background_button = self.findChild(QtWidgets.QRadioButton,"radioButton_ShowBackground")
		self.show_background_button.clicked.connect(self.graph_data_basic_plus_background)
		
		self.subtract_background_button = self.findChild(QtWidgets.QRadioButton,"radioButton_SubtractBackground")
		self.subtract_background_button.clicked.connect(self.graph_data_basic_subtract_background)
		return

	def get_background(self):
		self.background_type_1 = self.background_combo.currentText()
		self.background_type_2 = self.background_combo_2.currentText()
		self.background = self.calculate_background(self.data,self.background_range_1,self.background_range_2,self.background_type_1,self.background_type_2)
		self.add_background_to_graph()
		return
	
	def calculate_background(self,data,range_1,range_2,type_1,type_2):
		"""function calculates background based on the selected type"""
		background = xps_BackgroundWindow.set_background(data,range_1,type_1)
		if background is None:
			return None
		else:
			background_2 = xps_BackgroundWindow.set_background(data,range_2,type_2)
			background = xps_BackgroundWindow.combine_with_second(data,background,background_2,range_1,range_2,type_1,type_2)		
		return background
	
	def add_background_to_graph(self):
		"""function adds background to graph"""
		if self.show_background_button.isChecked() == True:
			self.graph_data_basic_plus_background()
		elif self.subtract_background_button.isChecked() == True:
			self.graph_data_basic_subtract_background()
		else:
			pass		
		return
	
	def set_background(data,indsx,background_type):
		"""function  sets background based on a selcted type"""
		#if statements that selects which type of background is calculated from the data
		background = np.array([np.nan]*len(data[:,1]))
		if background_type == "None":
			return None
		elif background_type == "Linear":
			background_temp = background_calculation.linear_background(data,indsx)
			background[indsx[0]:indsx[1]] = background_temp
		elif background_type == "Full Range Shirley":
			background = background_calculation._calculate_shirley_background_full_range(data[:,1])
		elif background_type == "Limited Range Shirley":
			background_temp = background_calculation._calculate_shirley_background_full_range(data[indsx[0]:indsx[1],1])
			background[indsx[0]:indsx[1]] = background_temp
		elif background_type == "tougaard":
			pass#background = background_calculation.tougaard_calculate(data[:,0],data[:,1])
		else:
			pass
		return background
	
	def combine_with_second(data,background,background_2,range_1,range_2,type_1,type_2):
		"""function  sets background based on a selcted type"""
		#if statements that selects which type of background is calculated from the data
		if type_2 == "None":
			return background
		else:
			#function to determine where background arrays overlap
			background_final = np.zeros(len(data[:,1]))
			background_final[range_1[0]:range_1[1]] = background[~np.isnan(background)]	
			background_final[range_2[0]:range_2[1]] = background_2[~np.isnan(background_2)]
			over_inds = xps_BackgroundWindow.determine_overlap(range_1,range_2)
			if over_inds is not None:
				#function to combine the two background arrays
				over_back = (background[over_inds] + background_2[over_inds])/2
				background_final[over_inds] = over_back
		return background_final
	
	def determine_overlap(range_1,range_2):
		"""function to determine where the ranges overlap"""
		if range_1[0] <= range_2[0]:
			min = range_2[0]
		else:
			min = range_1[0]
		if range_1[1] >= range_2[1]:
			max = range_2[1]
		else:
			max = range_1[1]
		if min > max:
			return  np.arange(max,min)
		elif max > min:
			return np.arange(min,max)
		else:
			return None
		
	def determine_background_1_indx(range_1,over_inds):
		"""function to determine where the ranges overlap"""
		print("range_1",range_1[1],over_inds)
		if range_1[0] < over_inds[0]:		
			back1ind1 = range_1[0]
			back1ind2 = over_inds[0]
		elif range_1[1] > over_inds[-1]+1:
			back1ind1 = over_inds[-1]+1
			back1ind2 = range_1[1]
		else:
			back1ind1,back1ind2 = 0,0
		return np.arange(back1ind1,back1ind2)
			
	def determine_background_2_indx(range_2,over_inds):
		"""function to determine where the ranges overlap"""
		if range_2[0] < over_inds[0]:
			back2ind1 = range_2[0]
			back2ind2 = over_inds[0]
		elif range_2[1] > over_inds[-1]+1:
			back2ind1 = over_inds[-1]+1
			back2ind2 = range_2[1]
		else:
			back2ind1,back2ind2 = 0,0		
		return np.arange(back2ind1,back2ind2)
	
	def select_range1(self):
		"""allows user to select range for first background calculation"""
		self.add_background_to_graph()
		self.span = SpanSelector(self.MplWidget.canvas.axes, self.onselect1, 'horizontal',
                                       useblit=True, interactive=True, drag_from_anywhere=True, props=dict(alpha=0.5, facecolor='red'))		
		return
	
	def select_range2(self):
		self.add_background_to_graph()
		"""allows user to select range for second background calculation"""
		self.span = SpanSelector(self.MplWidget.canvas.axes, self.onselect2, 'horizontal',
                                       useblit=True, interactive=True, drag_from_anywhere=True, props=dict(alpha=0.5, facecolor='blue'))
		return

	
	def onselect1(self, xmin, xmax):
		idx1,idx2 = (np.abs(self.data[:,0] - xmin)).argmin(),(np.abs(self.data[:,0] - xmax)).argmin()
		self.background_range_1[0],self.background_range_1[1] = idx1,idx2 
		return

	def onselect2(self, xmin, xmax):
		idx1,idx2 = (np.abs(self.data[:,0] - xmin)).argmin(),(np.abs(self.data[:,0] - xmax)).argmin()
		self.background_range_2[0],self.background_range_2[1] = idx1,idx2 
		return
	
	def graph_data_basic_plus_background(self):
		"""Function to plot chosen XPS plots data set"""
		#initalizes variables X & Y for plotting data
		X,Y = self.data[:,0],self.data[:,1]
		self.MplWidget.canvas.axes.clear()
		self.MplWidget.canvas.axes.plot(X[self.sind[0]:self.sind[1]],Y[self.sind[0]:self.sind[1]])
		if self.background is None or self.show_background_button.isChecked() == False:
			pass
		else:
			self.MplWidget.canvas.axes.plot(X[self.sind[0]:self.sind[1]],self.background[self.sind[0]:self.sind[1]])
			#draws
		self.MplWidget.canvas.draw()		
		return
	
	def graph_data_basic_subtract_background(self):
		"""Function to plot chosen XPS plots data set"""
		#initalizes variables X & Y for plotting data
		X,Y = self.data[:,0],self.data[:,1]
		self.MplWidget.canvas.axes.clear()
		if self.background is None or self.subtract_background_button.isChecked() == False:
			self.MplWidget.canvas.axes.plot(X[self.sind[0]:self.sind[1]],Y[self.sind[0]:self.sind[1]])
		else:
			self.MplWidget.canvas.axes.plot(X[self.sind[0]:self.sind[1]],Y[self.sind[0]:self.sind[1]]-self.background[self.sind[0]:self.sind[1]])
			#draws
		self.MplWidget.canvas.draw()		
		return

	def done1(self):
		self.close()
		return
	
class background_calculation:
	def linear_background(data,sind):
		"""Calculates a background by subtracting a linear line set between two points"""
		p1,p2 = np.array([data[sind[0],0],data[sind[0],1]]),np.array([data[sind[1]-1,0],data[sind[1]-1,1]])
		m = (p2[1]-p1[1])/(p2[0]-p1[0])
		b = p1[1]-m*p1[0]
		background = np.array([m*x+b for x in data[sind[0]:sind[1],0]])
		return background
	
	def _calculate_shirley_background_full_range(xps: np.ndarray, eps=1e-7, max_iters=50, n_samples=5) -> np.ndarray:
		"""Core routine for calculating a Shirley background on np.ndarray data. 
		TAKEN FROM PYARPES PACKAGE"""
		background = np.copy(xps)
		cumulative_xps = np.cumsum(xps, axis=0)
		total_xps = np.sum(xps, axis=0)

		rel_error = np.inf

		i_left = np.mean(xps[:n_samples], axis=0)
		i_right = np.mean(xps[-n_samples:], axis=0)

		iter_count = 0

		k = i_left - i_right
		for iter_count in range(max_iters):
			cumulative_background = np.cumsum(background, axis=0)
			total_background = np.sum(background, axis=0)

			new_bkg = np.copy(background)

			for i in range(len(new_bkg)):
				new_bkg[i] = i_right + k * (
					(total_xps - cumulative_xps[i] - (total_background - cumulative_background[i]))
					/ (total_xps - total_background + 1e-5)
				)

			rel_error = np.abs(np.sum(new_bkg, axis=0) - total_background) / (total_background)

			background = new_bkg

			if np.any(rel_error < eps):
				break

		if (iter_count + 1) == max_iters:
			warnings.warn(
				"Shirley background calculation did not converge "
				+ "after {} steps with relative error {}!".format(max_iters, rel_error)
			)
		return background

	