import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import SpanSelector
import pandas as pd

# Fixing random state for reproducibility
import numpy as np

def dat_handling(filename):
    with open(filename) as f:
        lines = f.readlines()
        text = "".join(lines)
    text = text.replace(" ","")
    
    jgap = np.array([])
    raw_data = np.array([])
    energy = np.array([])
    i0 = np.array([])
    
    for line in lines:
        if line[0].isdigit() == False:
            pass
        else:
            jgap = np.append(jgap,float(line.split()[0]))
            energy = np.append(energy,float(line.split()[1]))
            raw_data = np.append(raw_data,float(line.split()[2]))
            i0 = np.append(i0,float(line.split()[3]))
    data = np.column_stack((energy,raw_data*1000))
    return data

def dat_handling_I06(filename):
    with open(filename) as f:
        lines = f.readlines()
        text = "".join(lines)
    text = text.replace(" ","")
    energy = np.array([])
    XAS_pc = np.array([])  
    XAS_pc_std = np.array([])
    XAS_nc  = np.array([])
    XAS_nc_std  = np.array([])
    XAS_avg   = np.array([])
    XAS_avg_std  = np.array([])
    XMCD  = np.array([])
    XMCD_std = np.array([])
    
    for line in lines:
        if line[0].isdigit() == False:
            pass
        else:
            energy = np.append(energy,float(line.split()[0]))
            XAS_pc = np.append(XAS_pc,float(line.split()[1]))
            XAS_pc_std = np.append(XAS_pc_std,float(line.split()[2]))
            XAS_nc = np.append(XAS_nc,float(line.split()[3]))
            XAS_nc_std = np.append(XAS_nc_std,float(line.split()[4]))
            XAS_avg = np.append(XAS_avg,float(line.split()[5]))
            XAS_avg_std = np.append(XAS_avg_std,float(line.split()[6]))
            XMCD = np.append(XMCD,float(line.split()[7]))
            XMCD_std = np.append(XMCD_std,float(line.split()[8]))
            
    data1 = np.column_stack((energy,XAS_pc))
    data2 = np.column_stack((energy,XAS_nc))
    data3 = np.column_stack((energy,XAS_avg))
    data4 = np.column_stack((energy,XMCD))
    


    return data1,data2,data3,data4

def xy_handling(filename):
    with open(filename) as f:
        lines = f.readlines()
        text = "".join(lines)
    text = text.replace(" ","")
    x = np.array([])
    y = np.array([])
    for line in lines:
        if line[0].isdigit() == False:
            pass
        else:
            x = np.append(x,float(line.split()[0]))
            y = np.append(y,float(line.split()[1]))
    data = np.column_stack((x,y))
    return data
    
def dat_handling_IELTS(file):
    """function to open a .dat file """
    with open(file) as f:
        lines = f.readlines()
        text = "".join(lines)
    array = np.array([])
    i = 0
    for line in lines:
        i+=1
        if "[DATA]" in line:
            break
    df = pd.read_csv(file,dtype = float, skiprows=i+1,sep='\t',header=None)
    voltage = df[0]
    current_avg = df[1]
    didv_avg = df[2]
    di2dv2_avg = df[3]

    data = np.array([np.column_stack((voltage,di2dv2_avg))])
    return data




