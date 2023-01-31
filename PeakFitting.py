import lmfit as lm
import numpy as np
from xps_background import background_calculation
import matplotlib.pyplot as plt

class peak_fitting:
    """Peak fitting class provides functions for finding the best model to fit to the XPS data"""
    def peak_model(self,peak_no):
        """Function that creates an object array, each representing one peak"""
        #initializes an empty array, before populating each element with a VoightModel object
        #Each peak is labeled in the fit report by the prefix v{i}
        peaks = np.empty(peak_no,dtype = "object")
        for i in range(peak_no):
            peaks[i] = lm.models.VoigtModel(prefix = "v{}".format(i))
        return peaks

    def generate_model_guess(self,x,y,peaks,peak_no):
        pars = peaks[0].guess(y,x = x)
        #sets the variable model to contain the combine peak information
        model = peaks[0]
        for i in range(1,peak_no):
            pars.update(peaks[i].make_params())
            model = model + peaks[i]
        return model,pars

    def generate_model_with_params(self,x,y,peaks,peak_no,inital_peak_params):
        pars = peaks[0].make_params(center = x[inital_peak_params[0][0]])
        #sets the variable model to contain the combine peak information
        model = peaks[0]
        for i in range(1,peak_no):
            pars.update(peaks[i].make_params(center = x[inital_peak_params[0][i]]))
            model = model + peaks[i]
        return model,pars

    def contraints(self):
        return

    def fit_peaks(self,data,peak_no,inital_peak_params,background):
        """Function to fit the peak model to the XPS data"""
        #Initalize the data and the peak model into the variables x,y
        x,y = data[:,0],data[:,1]-background
        # initalizes the peak model into the variable peaks and guess some intial parameters to the first entry.
        peaks = self.peak_model(peak_no)
        if inital_peak_params is None:
            self.model,self.pars = self.generate_model_guess(x,y,peaks,peak_no)
        else:
            self.model,self.pars = self.generate_model_with_params(x,y,peaks,peak_no,inital_peak_params)
        #fits the peak model to the XPS data
        fit_results = self.model.fit(y,self.pars,x=x)
        print("\nFitting completed!\n.")
        return fit_results

    def plot_fitting(self,x,fit_results,background):
        """Function to display the fitting data, including the fit report following fitting calculation"""
        #plots the data and the fitted curve, along with printing the fit report
        self.MplWidget.canvas.axes.plot(x,fit_results.best_fit+background)
        self.MplWidget.canvas.draw()
        return

    def display_fit_report(self):
        print(self.fit_results.fit_report())
        return

    def plot_peak_model(self,x,fit_results,peak_no,background):
        comps = fit_results.eval_components()
        for i in range(peak_no):
            self.MplWidget.canvas.axes.fill(x,comps["v{}".format(i)]+background,alpha = 0.5)
        self.MplWidget.canvas.draw()
        return