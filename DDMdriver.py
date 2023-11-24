import numpy as np
import re
import pandas as pd
import matplotlib.pyplot as plt
import os
from scipy import signal
import subprocess
import time
import psutil
from JIN_pylib import Data2D_XT

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
        self.Reservoir_Center_Depth = 0
        self.Reservoir_Domain_X_Minimum = -200
        self.Reservoir_Domain_X_Maximum = 200
        self.Reservoir_Domain_Y_Minimum = -200
        self.Reservoir_Domain_Y_Maximum = 200
        self.Reservoir_Domain_Z_Minimum = -200
        self.Reservoir_Domain_Z_Maximum = 200

        # MONIT section
        self.well_path_X_begin = -200
        self.well_path_X_end = 200
        self.well_path_Y_begin = 25
        self.well_path_Y_end = 25
        self.well_path_Z_begin = 0
        self.well_path_Z_end = 0
        self.well_path_measured_depth_begin = 0
        self.gauge_length = 5.0

        self.output_at_timesteps = 35

        # data IO
        self.default_input_file = './input.am'
        self.default_output_file = './Outputs/MonitorWell_1.dat'
    
    def set_reservoir_Y_half_length(self, half_length):
        self.Reservoir_Domain_Y_Minimum = -half_length
        self.Reservoir_Domain_Y_Maximum = half_length
    
    def set_reservoir_X_half_length(self, half_length):
        self.Reservoir_Domain_X_Minimum = -half_length
        self.Reservoir_Domain_X_Maximum = half_length
    
    def set_reservoir_Z_half_length(self, half_length):
        self.Reservoir_Domain_Z_Minimum = -half_length
        self.Reservoir_Domain_Z_Maximum = half_length
    
    def set_well_distance(self, distance):
        self.well_path_Y_begin = distance
        self.well_path_Y_end = distance
    
    def set_well_depth(self, depth):
        self.well_path_Z_begin = depth
        self.well_path_Z_end = depth

    def set_slant_well(self, x0,y0,theta):
        a = np.tan(np.deg2rad(theta))
        b = y0-a*x0
        xmin = self.Reservoir_Domain_X_Minimum
        xmax = self.Reservoir_Domain_X_Maximum
        ymin = self.Reservoir_Domain_Y_Minimum
        ymax = self.Reservoir_Domain_Y_Maximum

        # calculate cross-points at boundaries:
        if a == 0:
            a = 1e-16
        pxs = np.array([xmin,xmax,(ymin-b)/a,(ymax-b)/a])
        pys = np.array([xmin*a+b, xmax*a+b,ymin,ymax])

        def in_range(x,y):
            if x>=xmin and x<=xmax and y>=ymin and y<=ymax:
                return True
            else:
                return False
        ind = [in_range(x,y) for x,y in zip(pxs,pys)]
        wellx=pxs[ind]
        welly=pys[ind]

        ind = np.argsort(wellx)
        wellx = wellx[ind]
        welly = welly[ind]

        self.well_path_X_begin = wellx[0]
        self.well_path_X_end = wellx[1]
        self.well_path_Y_begin = welly[0]
        self.well_path_Y_end = welly[1]

    
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

    def uniform_fracture(self,height, width):
        self.injection_time = _generate_numbers_string(self.start_time, self.end_time, self.number_of_time_steps)

        self.elements = _generate_elements(self.number_of_time_steps, self.number_of_fracture_branch, self.frac_element_half_length)

        width_all = []
        height_all = []
        xgrid_all = []

        for j in range(self.number_of_time_steps):
            current_length = (j+1)*self.frac_element_half_length*2
            elem_x = np.arange(j+1)*self.frac_element_half_length*2 + self.frac_element_half_length

            current_width = np.ones_like(elem_x)*width
            current_height = np.ones_like(elem_x)*height

            width_all.append(current_width)  
            height_all.append(current_height)
            xgrid_all.append(elem_x)
        
        self.xgrid_all = xgrid_all
        self.width_all = width_all
        self.height_all = height_all
    
    def run(self):
        self.output()
        run_model_process()
        check_model_run()
    
    def check_model_run(self):
        check_model_run()
    
    def get_result(self,dataset='Strain'):
        t,d,data = read_output(self.default_output_file,dataset=dataset)
        output_data = Data2D_XT.Data2D()
        output_data.taxis = t*60
        output_data.daxis = d
        output_data.data = data.T
        self.data = output_data
        return output_data
    
    def find_peaks(self):
        find_peaks_Data2D(self.data)
        return self.data.peaks_ind
    
    def plot_peaks(self,style='kx',**kwargs):
        self.data.plot_waterfall(**kwargs)
        plot_peaks_Data2D(self.data,style=style)
    
    def get_tip_dist(self, time_step):
        # calculate tip distance to each of the fracture element
        wellx = self.well_path_Y_begin
        welly = self.well_path_Z_begin
        dx = np.array(self.xgrid_all[time_step])-wellx
        dy = np.zeros_like(dx)
        for i in range(len(dy)):
            if welly < self.height_all[time_step][i]/2:
                dy[i] = 0
            else:
                dy[i] = welly - self.height_all[time_step][i]/2
        dists = np.sqrt(dx**2+dy**2)

        return np.min(dists)

    def output_peak_attributes(self):
        self.find_peaks()
        peak_ind = np.array(self.data.peaks_ind)
        results = []
        for i in range(self.data.data.shape[1]):
            tip_dist = self.get_tip_dist(i)
            if tip_dist < self.frac_element_half_length*2:
                break
            ind = peak_ind[:,0]==i
            peak_locs = self.data.daxis[peak_ind[ind,1]]
            peak_dist = np.max(peak_locs)-np.min(peak_locs)
            peak_amp = np.median(self.data.data[peak_ind[ind,1],peak_ind[ind,0]])
            frac_width = self.width_all[i][-1]
            frac_height = self.height_all[i][-1]
            well_depth = self.well_path_Z_begin
            results.append({'timestep':i,
                            'peak_dist':peak_dist,
                            'peak_amp':peak_amp,
                            'tip_dist':tip_dist,
                            'frac_width':frac_width,
                            'frac_height':frac_height,
                            'frac_elem_half_length':self.frac_element_half_length,
                            'well_depth':well_depth})
        results = pd.DataFrame(results)
        return results

            
            

            
    
    def draw_fracture(self,time_step_skip=10):
        plt.subplot(1,2,1)
        for i in range(0,self.number_of_time_steps,time_step_skip):
            x = self.xgrid_all[i]
            x = np.concatenate([x,[x[-1]]])
            y = self.height_all[i]/2
            y = np.concatenate(([y,[0]]))
            plt.plot(x,y)
        plt.xlabel('y')
        plt.ylabel('Half Height')

        plt.subplot(1,2,2)
        for i in range(0,self.number_of_time_steps,time_step_skip):
            x = self.xgrid_all[i]
            x = np.concatenate((x,[x[-1]]))
            y = self.width_all[i]
            y = np.concatenate(([y,[0]]))
            plt.plot(x,y,label=f'Time step: {i}')
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.xlabel('y')
        plt.ylabel('Width')


    
    def output(self,filename=None):
        if filename is None:
            filename = self.default_input_file

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
    numbers_str = "   ".join(["{:.1f}".format(num) for num in numbers])
    return numbers_str

def _generate_elements(number_of_time_steps, number_of_fracture_branch, half_length):
    radian_angle_positive = 1.570796
    radian_angle_negative = -1.570796
    x_mid = 0.000000

    entries = []
    for branch in range(1, number_of_fracture_branch + 1):
        y_mid = half_length
        if branch % 2 == 1:  # Odd branch
            radian_angle = radian_angle_positive
        else:  # Even branch
            radian_angle = radian_angle_negative
            y_mid = -y_mid

        for _ in range(number_of_time_steps):
            entries.append((x_mid, y_mid, half_length, radian_angle))
            y_mid += 2*half_length if radian_angle == radian_angle_positive else -2*half_length

    return entries


def read_output(filename, dataset='Strain'):
    """
    This function reads the output from a specified file.

    Parameters:
    filename (str): The name of the file to read.
    dataset (str, optional): The type of data to extract from the file. 
                             Options are  'U', 'Strain', 'Strain_Rate'
                             Default is 'Strain'.

    Returns: Time, Depth, data
    """
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
    s = np.asarray(d[dataset])
    data = s.reshape(len(Time)+1,len(Depth))
    data = data[1:,:] # remove the first row

    return Time, Depth, data

def check_model_run(input_file='./input.am',output_file='./Outputs/MonitorWell_1.dat'):
    if os.path.getmtime(input_file) > os.path.getmtime(output_file):
        raise Exception("input file is newer than result file, run exe first!")

def find_peaks_Data2D(DASdata):
    data = DASdata.data
    peaks = []
    skull = []
    for i in range (data.shape[1]):
        peak, _ = signal.find_peaks(data[:,i])
        peaks.append(peak)
        
    df = pd.DataFrame()
    df['loc'] = peaks

    for k in range(data.shape[1]):
        for j in range(df['loc'][k].shape[0]):
            skull.append([k,df['loc'][k][j]])
    
    DASdata.peaks_ind = skull
    
    return skull

def plot_peaks_Data2D(DASdata,style='kx'):
    x = []
    y = []
    for p in DASdata.peaks_ind:
        x.append(DASdata.taxis[p[0]])
        y.append(DASdata.daxis[p[1]])
    plt.plot(x,y,style)


def find_peaks(data):
    peaks = []
    skull = []
    for i in range (data.shape[0]):
        peak, _ = signal.find_peaks(data[i,:])
        peaks.append(peak)
        
    df = pd.DataFrame()
    df['loc'] = peaks

    for k in range(data.shape[0]):
        for j in range(df['loc'][k].shape[0]):
            skull.append([k,df['loc'][k][j]])
    
    return skull

def run_model_process():
    # Define the command and arguments
#    command = 'Fiber_Strain_rate_Engine.exe'
    command = '3D_DDM_slantedwells_04_15_22.exe'
    input_data = 'input.am'

    # Function to terminate a process and its children
    def terminate_process(proc):
        try:
            # Get the list of child processes
            child_procs = proc.children(recursive=True)

            # First, try to terminate the process
            proc.terminate()

            # Give it some time to terminate gracefully
            gone, alive = psutil.wait_procs([proc], timeout=3)

            # If the process is still alive, then kill it forcefully
            if alive:
                for p in alive:
                    p.kill()

            # Now, handle the child processes
            gone, alive = psutil.wait_procs(child_procs, timeout=3)
            if alive:
                for p in alive:
                    p.kill()

            return True
        except psutil.NoSuchProcess:
            return False  # The process does not exist
        except psutil.AccessDenied:
            return False  # Usually raised when requiring administrative privileges
        except Exception as e:
            print(f"An unexpected error occurred: {e}")
            return False

    # Record the start time
    start_time = time.time()

    # Start the process with suppressed output
    try:
        process = subprocess.Popen(
            f'echo {input_data} | {command}',
            shell=True,
            stdout=subprocess.DEVNULL,  # Suppress standard output
            stderr=subprocess.DEVNULL,  # Suppress standard error (optional, remove if error info is needed)
        )

        # Wait for a certain amount of time or until process finishes
        for _ in range(3):  # for example, check for 3 seconds
            if process.poll() is not None:
                break
            time.sleep(1)

        # If the process is still running, we need to terminate it and its children
        if process.poll() is None:
            proc = psutil.Process(process.pid)  # Get a psutil.Process instance for more features
            if terminate_process(proc):
                print("Process and its children have been terminated.")
            else:
                print("Failed to terminate the process (it may not exist, or require higher privileges).")
        else:
            print("Process had already ended.")

    except Exception as e:
        print(f"An error occurred: {e}")

    # Record the end time
    end_time = time.time()

    # Calculate and print the total execution time
    execution_time = end_time - start_time
    print(f"Total execution time: {execution_time:.2f} seconds")


def get_uniform_filename(fracture_half_height,fracture_width,well_depth):
    return 'results/uniform_H{}_W{}_WD{}.feather'.format(fracture_half_height,fracture_width,well_depth)

def get_ellipse_filename(fracture_half_height,fracture_width,well_depth):
    return 'results/ellipse_H{}_W{}_WD{}.feather'.format(fracture_half_height,fracture_width,well_depth)