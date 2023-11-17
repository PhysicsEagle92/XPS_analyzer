from nexusformat.nexus import *
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from contextlib import redirect_stdout


user = os.environ['USER']


def get_nexus_data(file):
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
    folder_path = f"/home/{user}/i09-jan-2023-data/si31574-1/"

    # INPUT: Detector entry ############################################
    # This is necessary because the nexus file can have multiple entries
    entry_string = "entry1"
    detector = "ew4000"
    file_list = [236656] # if plot_all_flag is True, this will be ignored, if True only these will be plotted.
    plot_all_flag = True


    # INPUT: Graph  ####################################################
    figure_size = (16, 10)  # size
    font_size = 36
    font_size_label = 20
    font_size_text = 20
    label_coefficient = 0.8
    marker_size = 30
    edge_width = 0.8
    marker_transparency = 0.15
    save_path = f"/home/{user}/i09-jan-2023-data/si31574-1/graphs/"

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
            data_list, metadata_list = get_nexus_data(file)
            plot_xps_data(data_list, metadata_list)
            #print(metadata_list)
            #plt.show()
        except:
            error_list.append(file_name)
            print(f"File {file_name} came out with an error")
error_path = f"{folder_path}/graphs"
write_list(error_list,error_path)



