from nexusformat.nexus import nxload
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import Slider, Button
import os

def get_file_list_from_folder(folder):
    """Function that returns a list of files in a folder"""
    file_list = os.listdir(folder)
    full_list = [os.path.join(folder,i) for i in file_list]
    return full_list

def check_polar(nxfile):
    """checks polarisation value"""
    polar = nxfile["entry"]["instrument"]["id"]["polarisation"].nxvalue
    return polar

def get_nexus_data_I06zacscan(file):
    """Function that loads the data from a nexus file and returns it as a list of numpy arrays"""
    entry_string = "entry"
    data_region_list = []
    metadata_region_list = None
    y_array = file[entry_string]["instrument"]["fesData"]["idio"].nxvalue
    x_array = file[entry_string]["instrument"]["fesData"]["iddenergy"].nxvalue
    data_region_list.append({"x": x_array, "y": y_array})
    return data_region_list, metadata_region_list

def get_nexus_data_I06stepscan(file):
    """Function that loads the data from a nexus file and returns it as a list of numpy arrays"""
    entry_string = "entry"
    data_region_list = []
    metadata_region_list = None
    y_array = file[entry_string]["ca61sr"]["ca61sr"].nxvalue
    x_array = file[entry_string]["ca61sr"]["energy"].nxvalue
    data_region_list.append({"x": x_array, "y": y_array})
    return data_region_list, metadata_region_list

def fix_max(ypos,y1):
    """Function that fixes the maximum of two arrays"""
    y1_max = np.max(y1)
    pos1 = np.where(y1 == y1_max)[0][0]
    if ypos > pos1:
        y1 = np.roll(y1,ypos-pos1)
    elif ypos < pos1:
        y1 = np.roll(y1,ypos-pos1)
    return y1

def combine_spectra(spectra_list):
    spec = nxload(spectra_list[0])
    data = get_nexus_data_I06zacscan(spec)
    spectra_length = len(data[0][0]['x'])
    #ypos = np.where(data[0][0]['y'] == np.max(data[0][0]['y']))[0][0]
    x = data[0][0]['x']
    pc,nc = np.zeros(spectra_length),np.zeros(spectra_length)
    pc_count,nc_count = 0,0
    for spectra in spectra_list:
        spec = nxload(spectra)
        data = get_nexus_data_I06zacscan(spec)
        #data[0][0]['y'] = fix_max(ypos,data[0][0]['y'])
        if check_polar(spec) == "pc":
            pc = pc + data[0][0]['y']
            pc_count +=1
        else:
            nc = nc + data[0][0]['y']
            nc_count+=1
    return x,pc/pc_count,nc/nc_count

def graph_slider(x1,x2,y1,y2):
    """function that draws two plots with sliders that shift each plot in x direction"""
    fig, (ax1, ax2) = plt.subplots(2,1)
    line1, = ax1.plot(x1,y1)
    line2, = ax1.plot(x2,y2)
    line3, = ax2.plot(x1,y1-y2)
    #add slider to graph
    axcolor = 'lightgoldenrodyellow'
    axpos = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
    spos = Slider(axpos,'slide', -10, 10, valinit=0, valstep=1)
    
    def update(val):
        pos = spos.val
        line1.set_ydata(y1)
        line2.set_ydata(np.roll(y2,int(pos)))
        
        line3.set_ydata(y1-np.roll(y2,int(pos)))
        fig.canvas.draw_idle()
        
    spos.on_changed(update)
    plt.show()
    return
    
    

folder = "C:\\Users\ppxcf1\\Desktop\\Data\\I06-beamtime-July\\Down_Set"
file_list = get_file_list_from_folder(folder)
x,pc,nc = combine_spectra(file_list)
graph_slider(x,x,pc,nc)





