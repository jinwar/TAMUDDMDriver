import numpy as np
import re
import pandas as pd
import matplotlib.pyplot as plt
import os


class DDM:


    def __init__(self):
        self.default_values()
        pass;


    def default_values(self):
        # stimulation
        self.frac_element_half_length = 1
        self.number_of_time_steps = 60
        self.start_time = 0.0  # Example value, can be changed
        self.end_time = 59.0   # Example value, can be changed

        # Constants and static values
        self.number_of_fracture_branch = 2
        self.number_of_element_in_Z_direction = 1
        self.number_of_dimension = 3
        self.number_of_monitor_wells = 1
        self.number_of_element_each_frac_branch =  ", ".join(map(str, tuple([self.number_of_time_steps] * self.number_of_fracture_branch)))
        self.number_of_monitor_point_each_well = 200
        self.Young_Modulus = 21400000000.0
        self.Poisson_Ratio = 0.26
        
        # RESVR section
        self.Reservoir_Center_Depth = 2000
        self.Reservoir_Domain_X_Minimum = -150
        self.Reservoir_Domain_X_Maximum = 150
        self.Reservoir_Domain_Y_Minimum = -150
        self.Reservoir_Domain_Y_Maximum = 150
        self.Reservoir_Domain_Z_Minimum = 1900
        self.Reservoir_Domain_Z_Maximum = 2100

        # MONIT section
        self.well_path_X_begin = -150
        self.well_path_X_end = 150
        self.well_path_Y_begin = 25
        self.well_path_Y_end = 25
        self.well_path_Z_begin = 2000
        self.well_path_Z_end = 2000
        self.well_path_measured_depth_begin = 0
        self.gauge_length = 5.0

        self.output_at_timesteps = 35

    
    def ellipse_fracture(self,h_ratio, w_ratio):
        self.injection_time = _generate_numbers_string(self.start_time, self.end_time, self.number_of_time_steps)

        self.elements = _generate_elements(self.number_of_time_steps, self.number_of_fracture_branch, self.frac_element_half_length)

        width_all = []
        height_all = []
        xgrid_all = []

        for j in range(self.number_of_time_steps):
            current_length = (j+1)*self.frac_element_half_length*2
            elem_x = np.arange(j+1)*self.frac_element_half_length*2 + self.frac_element_half_length

            current_width = np.sqrt((current_length**2-elem_x**2)*w_ratio**2)
            current_height = np.sqrt((current_length**2-elem_x**2)*h_ratio**2)

            width_all.append(current_width)  
            height_all.append(current_height)
            xgrid_all.append(elem_x)
        
        self.xgrid_all = xgrid_all
        self.width_all = width_all
        self.height_all = height_all
    
    def draw_fracture(self,time_step_skip=10):
        plt.subplot(1,2,1)
        for i in range(0,self.number_of_time_steps,time_step_skip):
            plt.plot(self.xgrid_all[i],self.height_all[i])
        plt.xlabel('y')
        plt.ylabel('Height')

        plt.subplot(1,2,2)
        for i in range(0,self.number_of_time_steps,time_step_skip):
            plt.plot(self.xgrid_all[i],self.width_all[i])
        plt.xlabel('y')
        plt.ylabel('Width')


    
    def output(self,filename):

        with open(filename, 'w') as file:
            # Writing the MEMORY section
            file.write("Input file for strain calculation\n")
            file.write(">>>MEMORY\n")
            file.write("&Basic_Parameter_Definitions           number_of_fracture_branch            = {:>5}\n".format(self.number_of_fracture_branch))
            file.write("                                       number_of_element_in_Z_direction     = {:>5}\n".format(self.number_of_element_in_Z_direction))
            file.write("                                       number_of_dimension                  = {:>5}\n".format(self.number_of_dimension))
            file.write("                                       number_of_time_steps                 = {:>5}\n".format(self.number_of_time_steps))
            file.write("                                       number_of_monitor_wells              = {:>5}\n".format(self.number_of_monitor_wells))
            file.write("       /\n")
            file.write("&Number_Attribute_Per_Segment          number_of_element_each_frac_branch   = {}\n".format(self.number_of_element_each_frac_branch))
            file.write("\n")
            file.write("                                       number_of_monitor_point_each_well    = {:>5}\n".format(self.number_of_monitor_point_each_well))
            file.write("       /\n")
            file.write("<<<\n")
            
            # Writing the ROCKS section
            file.write(">>>ROCKS\n")
            file.write("&Rock_Mechanical_Properties           Young_Modulus   = {:>14}\n".format(self.Young_Modulus))
            file.write("                                      Poisson_Ratio    = {:>14}\n".format(self.Poisson_Ratio))
            file.write("       /\n")
            file.write("<<<\n")
            
            # Writing the TIMES section
            file.write(">>>TIMES\n")
            file.write("&Injecton_Times                       injection_time   = {}\n".format(self.injection_time))
            file.write("       /\n")
            file.write("<<<\n")
            
            # Writing the RESVR section
            file.write(">>>RESVR\n")
            file.write("&Reservoir_Condition_Parameters     Reservoir_Center_Depth                  = {:>14}\n".format(self.Reservoir_Center_Depth))
            file.write("                                      Reservoir_Domain_X_Minimum             = {:>14}\n".format(self.Reservoir_Domain_X_Minimum))
            file.write("                                      Reservoir_Domain_X_Maximum             = {:>14}\n".format(self.Reservoir_Domain_X_Maximum))
            file.write("                                      Reservoir_Domain_Y_Minimum             = {:>14}\n".format(self.Reservoir_Domain_Y_Minimum))
            file.write("                                      Reservoir_Domain_Y_Maximum             = {:>14}\n".format(self.Reservoir_Domain_Y_Maximum))
            file.write("                                      Reservoir_Domain_Z_Minimum             = {:>14}\n".format(self.Reservoir_Domain_Z_Minimum))
            file.write("                                      Reservoir_Domain_Z_Maximum             = {:>14}\n".format(self.Reservoir_Domain_Z_Maximum))        
            file.write("       /\n")
            file.write("<<<\n")
            
            # Writing the MONIT section
            file.write(">>>MONIT\n")
            file.write("&Monitor_Wells_Definitions          well_path_X_begin                           = {:>14}\n".format(self.well_path_X_begin))
            file.write("                                      well_path_X_end                           = {:>14}\n".format(self.well_path_X_end))
            file.write("                                      well_path_Y_begin                         = {:>14}\n".format(self.well_path_Y_begin))
            file.write("                                      well_path_Y_end                           = {:>14}\n".format(self.well_path_Y_end))
            file.write("                                      well_path_Z_begin                         = {:>14}\n".format(self.well_path_Z_begin))
            file.write("                                      well_path_Z_end                           = {:>14}\n".format(self.well_path_Z_end))
            file.write("                                      well_path_measured_depth_begin            = {:>14}\n".format(self.well_path_measured_depth_begin))
            file.write("                                      gauge_length                              = {:>14}\n".format(self.gauge_length))
            file.write("       /\n")
            file.write("<<<\n")
            
            # Writing the OUTPUT section
            file.write(">>>OUTPUT\n")
            file.write("&Specified_Output_Timesteps          output_at_timesteps                        = {:>14}\n".format(self.output_at_timesteps))
            file.write("       /\n")
            file.write("<<<\n")
            
            # Writing the ELEME section
            file.write(">>>ELEME-----X_mid,       Y_mid,          half_length,    radian_angle\n")
            for elem in self.elements:
                file.write("       {:>14.6f}   {:>14.6f}   {:>14.6f}   {:>14.6f}\n".format(elem[0], elem[1], elem[2], elem[3]))

            # Writing displacement section
            file.write("::: Displacement Discontinuity of Each Element at Each Time Step\n")
            for i in range(self.number_of_time_steps):
                for j in range(len(self.width_all[i])):
                    line_values = [self.width_all[i][j], 0, self.Reservoir_Center_Depth - self.height_all[i][j]/2, self.Reservoir_Center_Depth + self.height_all[i][j]/2]
                    if self.number_of_fracture_branch % 2 == 0:
                        line_values = line_values + line_values  # double the columns for even branches
                        file.write("   ".join(["{:>14.6f}".format(val) for val in line_values]) + "\n")

            # Closing section and end of file
            file.write("<<<\n")
            file.write("ENDCY----1----*----2----*----3----*----4----*----5----*----6----*----7----*----8\n")

        print(f"File saved at: {filename}")
    

# define help functions

def _generate_numbers_string(start, end, num_points):
    numbers = np.linspace(start, end, num_points)
    numbers_str = "   ".join(map(str, numbers))
    return numbers_str

def _generate_elements(number_of_time_steps, number_of_fracture_branch, half_length):
    radian_angle_positive = 1.570796
    radian_angle_negative = -1.570796
    x_mid = 0.000000

    entries = []
    for branch in range(1, number_of_fracture_branch + 1):
        y_mid = 0.5
        if branch % 2 == 1:  # Odd branch
            radian_angle = radian_angle_positive
        else:  # Even branch
            radian_angle = radian_angle_negative
            y_mid = -y_mid

        for _ in range(number_of_time_steps):
            entries.append((x_mid, y_mid, half_length, radian_angle))
            y_mid += 2*half_length if radian_angle == radian_angle_positive else -2*half_length

    return entries


def read_output(filename):
    # Read the content of the file
    with open(filename, 'r') as file:
        content = file.read()

    # Use a regex pattern to find the value of N
    match = re.search(r'ZONE N =\s*(\d+)', content)

    # Extract the value if found
    if match:
        N_value = int(match.group(1))
        print(f"N = {N_value}")
    else:
        print("N value not found in the file.")

    data = np.loadtxt(filename, skiprows = 3,max_rows = N_value, unpack = True)

    d = pd.DataFrame(data.T)
    d.columns = ['Time', 'Depth', 'U', 'Strain', 'Strain_Rate']

    Time = d['Time'].unique()
    Depth = d['Depth'].unique()
    s = np.asarray(d['Strain'])
    data = s.reshape(len(Time)+1,len(Depth))

    return Time, Depth, data

def check_model_run(input_file,output_file):
    if os.path.getmtime(input_file) > os.path.getmtime(output_file):
        raise Exception("input file is newer than result file, run exe first!")