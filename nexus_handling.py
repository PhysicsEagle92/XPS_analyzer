from nexusformat.nexus import nxload
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import os
from contextlib import redirect_stdout


#user = os.environ['USER']

def get_nexus_data(file):
    """Function that loads the data from a nexus file and returns it as a list of numpy arrays"""
    entry_string = "entry1"
    if check_Spectra_vs_XAS(file,entry_string) == True:
        data_region_list, metadata_region_list = get_nexus_data_spectra(file)
    elif check_XAS_vs_ResPes(file,entry_string) == True:       
        data_region_list, metadata_region_list = get_nexus_data_ResPes(file)
    elif check_I06(file,entry_string):
        data_region_list, metadata_region_list = get_nexus_data_I06(file)
    else:
        data_region_list, metadata_region_list = get_nexus_data_XAS(file)
    return data_region_list, metadata_region_list

def check_Spectra_vs_XAS(file,entry_string):
    """Function that checks if the nexus file is a spectra file or an XAS file and returns the appropriate data"""
    if "ew4000" in file[entry_string]["instrument"]:
        return True
    else:
        return False
    
def check_XAS_vs_ResPes(file,entry_string):
    """Function that checks if the nexus file is a XAS file or a Resonant PES file and returns the appropriate data"""
    if "Res_Auger_N1s" in file[entry_string]:
        return True
    else:
        return False

def check_I06(file,entry_string):
    """Function that checks if the nexus file is a I06 file and returns the appropriate data"""
    if "fesData" in file[entry_string]:
        return True
    else:
        return False

def get_nexus_data_I06(file):
    """Function that loads the data from a nexus file and returns it as a list of numpy arrays"""
    entry_string = "entry"
    data_region_list = []
    metadata_region_list = None
    x_array = file[entry_string]["instrument"]["fastEnergy"]["value"].nxvalue
    y_array = file[entry_string]["instrument"]["fesData"]["C1"].nxvalue
    data_region_list.append({"x": x_array, "y": y_array})
    return data_region_list, metadata_region_list

def get_nexus_data_ResPes(file):
    """Function that loads the data from a nexus file and returns it as a list of numpy arrays"""
    entry_string = "entry1"
    detector = "Res_Auger_N1s"
    start_time = file[entry_string]['start_time'].nxvalue
    end_time = file[entry_string]['end_time'].nxvalue
    i0 = file[entry_string]["instrument"]["smpm"]["smpmiamp39"].nxvalue
    metadata_region_list = []
    data_region_list = []
    y_array = file[entry_string][detector]["spectrum_data"].nxvalue
    x_array = file[entry_string]["instrument"]['jenergy']["jenergy"].nxvalue*1000
    data_region_list.append({"x": x_array, "y": y_array})
    metadata_region_list.append({"i0": i0, "start_time": start_time, "end_time": end_time})
    return data_region_list, metadata_region_list


def get_nexus_data_XAS(file):
    """Function that loads the data from a nexus file and returns it as a list of numpy arrays"""
    entry_string = "entry1"
    detector = "scaler2"
    start_time = file[entry_string]['start_time'].nxvalue
    end_time = file[entry_string]['end_time'].nxvalue
    i0 = file[entry_string]["scaler2"]["smpmamp39"].nxvalue
    i01 = file[entry_string]["scaler2"]["sm5amp8"].nxvalue
    metadata_region_list = []
    data_region_list = []   
    y_array = file[entry_string]["instrument"][detector]["smpmamp39"].nxvalue
    x_array = file[entry_string]["instrument"]['jenergy']["jenergy"].nxvalue*1000
    data_region_list.append({"x": x_array, "y": y_array})
    metadata_region_list.append({"i0" : i0,"i01" : i01,"start_time" : start_time, "end_time" : end_time})
    data_region_list[0]["y"] = (data_region_list[0]["y"] / metadata_region_list[0]["i01"])*1000
    return data_region_list, metadata_region_list

def get_nexus_data_spectra(file):
    entry_string = "entry1"
    detector = "ew4000"
    start_time = file[entry_string]['start_time'].nxvalue
    end_time = file[entry_string]['end_time'].nxvalue
    region_name_list = file[entry_string]["instrument"][detector]["region_list"].nxvalue
    region_name_list = region_name_list.split(",")
    metadata_region_list = []
    data_region_list = []



    for region in region_name_list:
        attributes = file[entry_string]["instrument"][region]
        y_array = file[entry_string]["instrument"][region].spectrum_data.nxvalue  # Y data
        x_array = file[entry_string]["instrument"][region].energies.nxvalue  # X data
        y_array = y_array[0]
        x_array = x_array[0]

        acquisition_mode = file[entry_string]["instrument"][region].acquisition_mode.nxvalue
        angles = file[entry_string]["instrument"][region].angles.nxvalue
        energy_mode = file[entry_string]["instrument"][region].energy_mode.nxvalue
        energy_step = file[entry_string]["instrument"][region].energy_step.nxvalue
        excitation_energy = file[entry_string]["instrument"][region].excitation_energy.nxvalue
        external_io_data = file[entry_string]["instrument"][region].external_io_data.nxvalue
        fixed_energy = file[entry_string]["instrument"][region].fixed_energy.nxvalue
        high_energy = file[entry_string]["instrument"][region].high_energy.nxvalue
        image_data = file[entry_string]["instrument"][region].image_data.nxvalue
        lens_mode = file[entry_string]["instrument"][region].lens_mode.nxvalue
        local_name = file[entry_string]["instrument"][region].local_name.nxvalue
        low_energy = file[entry_string]["instrument"][region].low_energy.nxvalue
        number_of_iterations = file[entry_string]["instrument"][region].number_of_iterations.nxvalue
        number_of_slices = file[entry_string]["instrument"][region].number_of_slices.nxvalue
        pass_energy = file[entry_string]["instrument"][region].pass_energy.nxvalue
        step_time = file[entry_string]["instrument"][region].step_time.nxvalue
        total_steps = file[entry_string]["instrument"][region].total_steps.nxvalue
        total_time = file[entry_string]["instrument"][region].total_time.nxvalue

        #print(f"Metadata found are for file {file_name}: acquisition mode: {acquisition_mode}, angle: {angles}, energy mode: {energy_mode}, energy step: {energy_step}, excitation energy: {excitation_energy}, fixed energy: {fixed_energy}, high energy: {high_energy}, lens mode: {lens_mode}, local name: {local_name}, low energy: {low_energy}, number of iterations: {number_of_iterations}, number of slices: {number_of_slices}, pass energy: {pass_energy}, step time: {step_time}, total steps: {total_steps}, total time: {total_time}")

        metadata_region_list.append({"acquisition_mode": acquisition_mode, "angles": angles, "energy_mode": energy_mode,
                    "energy_step": energy_step, "excitation_energy": excitation_energy, "fixed_energy": fixed_energy,
                    "high_energy": high_energy, "lens_mode": lens_mode, "local_name": local_name,
                    "low_energy": low_energy, "number_of_iterations": number_of_iterations,
                    "number_of_slices": number_of_slices, "pass_energy": pass_energy, "step_time": step_time,
                    "total_steps": total_steps, "total_time": total_time, "start_time": start_time, "end_time": end_time, "x_array": x_array,
                    "y_array": y_array, "region_name": region, "attributes": attributes})
        data_region_list.append({"x": x_array, "y": y_array, "image_data": image_data, "external_io_data": external_io_data})

    return data_region_list, metadata_region_list


def plot_xps_data(data_list, metadata_list):

    for i, _ in enumerate(data_list):
        x = data_list[i]["x"]
        y = data_list[i]["y"]
        region = metadata_list[i]["region_name"]
        acquisition_mode = metadata_list[i]["acquisition_mode"]
        iterations = metadata_list[i]["number_of_iterations"]
        excitation_energy = metadata_list[i]["excitation_energy"]
        pass_energy = metadata_list[i]["pass_energy"]
        start_time = metadata_list[i]["start_time"].split("T")
        date = start_time[0]
        time = start_time[1].split(".")[0]

        figure_title = f"XPS Spectrum for {file_name}"

        fig = plt.subplots(figsize=figure_size)
        plt.title(figure_title, fontsize=font_size)
        plt.xlabel("Energy (eV)", fontsize=font_size_label)
        plt.ylabel("Counts [Arb. Units]", fontsize=font_size_label)
        plt.xticks(fontsize=font_size_label)
        plt.yticks(fontsize=font_size_label)
        sns.lineplot(x=x, y=y, color="#845EC2")
        sns.scatterplot(x=x, y=y, s=marker_size, edgecolor="#845EC2", facecolor="None", linewidth=edge_width, alpha=marker_transparency)
        text_string = (f"Region: {region} \n"
                       f"Excitation Energy: {round(excitation_energy, 2)} eV\n"
                       f"Pass Energy: {round(pass_energy, 2)} eV\n"
                       f"Acquisition Mode: {acquisition_mode}\n"
                       f"Iterations: {iterations}\n"
                       f"Date: {date}\n"
                       f"Time: {time}")
        plt.annotate(text_string, xy=(0.65, 0.65), xycoords="axes fraction", fontsize=font_size_text, color="#808080")
        figure_name = f"{file_name.split('.nxs')[0]}_{region}.png"
        plt.savefig(f"{save_path}{figure_name}", dpi=300, bbox_inches="tight")
        plt.close()
        print(f"Figure saved as {figure_name}")

def write_list(list,folder_path):
    n_list = ["{}\n".format(i) for i in list]
    with open(f'{folder_path}/files_with_errors.txt', 'w') as fp:
        fp.writelines(n_list)



if __name__ == "__main__":

    # INPUT: path to the nexus fil e############################################
    prefix = f"i09-"
    folder_path = f"C:/Users/ppxcf1/Desktop/Data/I09-beamtime-January/si31574-1/"

    # INPUT: Detector entry ############################################
    # This is necessary because the nexus file can have multiple entries
    entry_string = "entry1"
    detector = "ew4000"
    file_list = [236656] # if plot_all_flag is True, this will be ignored, if True only these will be plotted.
    plot_all_flag = False


    # INPUT: Graph  ####################################################
    figure_size = (16, 10)  # size
    font_size = 36
    font_size_label = 20
    font_size_text = 20
    label_coefficient = 0.8
    marker_size = 30
    edge_width = 0.8
    marker_transparency = 0.15
    save_path ="C:/Users/ppxcf1/Desktop/Data/I09-beamtime-January/saved_data"

    ### Check if save_path exists, if not create it
    if not os.path.exists(save_path):
        os.makedirs(save_path)

    # Data and Metadata #########################################################

    if plot_all_flag is True:
        file_list = []
        print(folder_path)
        for _fstring in os.listdir(folder_path):
            if _fstring.endswith(".nxs"):
                file_list.append(_fstring.strip(".nxs")[-6:])
    print(file_list)
    print(len(file_list))
    print(f"Files to be processed: {file_list}")

    error_list = []

    for id in file_list:
        file_name = f"{prefix}{id}.nxs"
        full_path = f"{folder_path}{file_name}"

        try:
            file = nxload(full_path)
            data_list, metadata_list = get_nexus_data(file,entry_string)
            plot_xps_data(data_list, metadata_list)
            #print(metadata_list)
            #plt.show()
        except:
            error_list.append(file_name)
            print(f"File {file_name} came out with an error")
#error_path = f"{folder_path}/graphs"
#write_list(error_list,error_path)

from nexusformat.nexus import nxload
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def plot_xps_data(data_list, metadata_list):

    for i, _ in enumerate(data_list):
        x = data_list[i]["x"]
        y = data_list[i]["y"]
        region = metadata_list[i]["region_name"]
        acquisition_mode = metadata_list[i]["acquisition_mode"]
        iterations = metadata_list[i]["number_of_iterations"]
        excitation_energy = metadata_list[i]["excitation_energy"]
        pass_energy = metadata_list[i]["pass_energy"]
        start_time = metadata_list[i]["start_time"].split("T")
        date = start_time[0]
        time = start_time[1].split(".")[0]

        figure_title = f"XPS Spectrum for {file_name}"

        fig = plt.subplots(figsize=figure_size)
        plt.title(figure_title, fontsize=font_size)
        plt.xlabel("Energy (eV)", fontsize=font_size_label)
        plt.ylabel("Counts [Arb. Units]", fontsize=font_size_label)
        plt.xticks(fontsize=font_size_label)
        plt.yticks(fontsize=font_size_label)
        sns.lineplot(x=x, y=y, color="#845EC2")
        sns.scatterplot(x=x, y=y, s=marker_size, edgecolor="#845EC2", facecolor="None", linewidth=edge_width, alpha=marker_transparency)
        text_string = (f"Region: {region} \n"
                       f"Excitation Energy: {round(excitation_energy, 2)} eV\n"
                       f"Pass Energy: {round(pass_energy, 2)} eV\n"
                       f"Acquisition Mode: {acquisition_mode}\n"
                       f"Iterations: {iterations}\n"
                       f"Date: {date}\n"
                       f"Time: {time}")
        plt.annotate(text_string, xy=(0.65, 0.65), xycoords="axes fraction", fontsize=font_size_text, color="#808080")
        figure_name = f"{file_name.split('.nxs')[0]}_{region}.png"
        plt.savefig(f"{save_path}{figure_name}", dpi=300, bbox_inches="tight")
        plt.close()
        print(f"Figure saved as {figure_name}")


def get_nexus_data2(file, entry_string="entry1", detector="ew4000"):
    try:
        start_time = file[entry_string]["start_time"].nxvalue
    except:
        print("Could not get the start time")
        start_time = ""
    try:
        end_time = file[entry_string]['end_time'].nxvalue
    except:
        print("Could not get the end time")
        end_time = ""
    try:
        region_name_list = file[entry_string]["instrument"][detector]["region_list"].nxvalue
        region_name_list = region_name_list.split(",")
    except:
        print("Could not get the region list")
        region_name_list = []
        # Get the region which is the first parameter after entry1 in the Nexus file structure
        temp_region = list(file[entry_string].keys())[0]
        region_name_list.append(temp_region)

    metadata_region_list = []
    data_region_list = []

    for region in region_name_list:
        try:
            attributes = file[entry_string]["instrument"][region]
        except:
            print("Could not get the attributes")
            attributes = ""
        try:
            spectrum_data = file[entry_string]["instrument"][region].spectrum_data.nxvalue  # Y data
        except:
            print("Could not get the spectrum data")
            spectrum_data = ""
        try:
            energies = file[entry_string]["instrument"][region].energies.nxvalue  # X data
        except:
            print("Could not get the energies")
            energies = ""
        try:
            acquisition_mode = file[entry_string]["instrument"][region].acquisition_mode.nxvalue
        except:
            print("Could not get the acquisition mode")
            acquisition_mode = ""
        try:
            angles = file[entry_string]["instrument"][region].angles.nxvalue
        except:
            angles = ""
        try:
            energy_mode = file[entry_string]["instrument"][region].energy_mode.nxvalue
        except:
            energy_mode = ""
        try:
            energy_step = file[entry_string]["instrument"][region].energy_step.nxvalue
        except:
            energy_step = ""
        try:
            excitation_energy = file[entry_string]["instrument"][region].excitation_energy.nxvalue
        except:
            excitation_energy = ""
        try:
            fixed_energy = file[entry_string]["instrument"][region].fixed_energy.nxvalue
        except:
            fixed_energy = ""
        try:
            high_energy = file[entry_string]["instrument"][region].high_energy.nxvalue
        except:
            high_energy = ""
        try:
            image_data = file[entry_string]["instrument"][region].image_data.nxvalue
        except:
            image_data = ""
        try:
            lens_mode = file[entry_string]["instrument"][region].lens_mode.nxvalue
        except:
            lens_mode = ""
        try:
            local_name = file[entry_string]["instrument"][region].local_name.nxvalue
        except:
            local_name = ""
        try:
            low_energy = file[entry_string]["instrument"][region].low_energy.nxvalue
        except:
            low_energy = ""
        try:
            number_of_iterations = file[entry_string]["instrument"][region].number_of_iterations.nxvalue
        except:
            number_of_iterations = ""
        try:
            number_of_slices = file[entry_string]["instrument"][region].number_of_slices.nxvalue
        except:
            number_of_slices = ""
        try:
            pass_energy = file[entry_string]["instrument"][region].pass_energy.nxvalue
        except:
            pass_energy = ""
        try:
            step_time = file[entry_string]["instrument"][region].step_time.nxvalue
        except:
            step_time = ""
        try:
            total_steps = file[entry_string]["instrument"][region].total_steps.nxvalue
        except:
            total_steps = ""
        try:
            total_time = file[entry_string]["instrument"][region].total_time.nxvalue
        except:
            total_time = ""
        try:
            external_io_data = file[entry_string]["scaler2"]["sm5amp8"].nxvalue
        except:
            print("Could not get the external_io_data")
            external_io_data = ""

        #print(f"Metadata found are for file {file_name}: acquisition mode: {acquisition_mode}, angle: {angles}, energy mode: {energy_mode}, energy step: {energy_step}, excitation energy: {excitation_energy}, fixed energy: {fixed_energy}, high energy: {high_energy}, lens mode: {lens_mode}, local name: {local_name}, low energy: {low_energy}, number of iterations: {number_of_iterations}, number of slices: {number_of_slices}, pass energy: {pass_energy}, step time: {step_time}, total steps: {total_steps}, total time: {total_time}")

        metadata_region_list.append({"acquisition_mode": acquisition_mode, "angles": angles, "energy_mode": energy_mode,
                    "energy_step": energy_step, "excitation_energy": excitation_energy, "fixed_energy": fixed_energy,
                    "high_energy": high_energy, "lens_mode": lens_mode, "local_name": local_name,
                    "low_energy": low_energy, "number_of_iterations": number_of_iterations,
                    "number_of_slices": number_of_slices, "pass_energy": pass_energy, "step_time": step_time,
                    "total_steps": total_steps, "total_time": total_time, "start_time": start_time, "end_time": end_time, "energies": energies,
                    "spectrum_data": spectrum_data, "region_name": region, "attributes": attributes})
        data_region_list.append({"energies": energies, "spectrum_data": spectrum_data, "image_data": image_data, "i0": external_io_data})

    return data_region_list, metadata_region_list
