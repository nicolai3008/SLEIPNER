# some_file.py
import sys
import numpy as np
import matplotlib.pyplot as plt
import os
from subprocess import call
# caution: path[0] is reserved for script path (or '' in REPL)
sys.path.insert(1, '../../')
from functions import *
import time
import shutil



def Theta2(L0,d):
    return 2*np.arcsin(L0/(2*d))
    


def read_previous(n = -1):
    with open(f'path_{order}.csv', 'r') as f:
        lines = f.readlines()
        previous = lines[n]
    previous = previous.split(';')
    return previous
   
   
def write_path(next):
    with open(f'path_{order}.csv', 'a') as f:
        f.write("\n")
        line = ""
        for key in next.keys():
            if key == "Variable":
                continue
            if type(next[key]) == list:
                line += ",".join([str(i) for i in next[key]]) + ";"
            else:
                line += str(next[key]) + ";"
        line += next["Variable"]+";"
        f.write(line)
    return True     

def parameters():
    # read previous
    previous = read_previous()
    # dictionary of parameters
    par_names = ["run", "L1", "L2", "L3", "mos1", "mos2", "col_h", "col_v"]
    d = {}
    i = 0
    for name, prev in zip(par_names,previous):
        if len(prev.split(',')) > 1:
            d[name] = [float(i) for i in prev.split(',')]
            d["Variable"] = name
        else:
            if i < 1:
                d[name] = int(prev)
            else:
                d[name] = float(prev)
        i += 1
    return d

def prev_parameters(var):
    # Read 2nd previous
    previous = read_previous(-2)
    # get previous value of variable
    par_names = ["run", "L1", "L2", "L3", "mos1", "mos2", "col_h", "col_v"]
    d = {}
    i = 0
    for name, prev in zip(par_names,previous):
        if len(prev.split(',')) > 1:
            d[name] = [float(i) for i in prev.split(',')]
        else:
            if i < 1:
                d[name] = int(prev)
            else:
                d[name] = float(prev)
        i += 1
    return d[var]
    

def read_data(folder, d):
    er = 0
    # read data
    run = d["run"]
    var = d["Variable"]
    intensities = []
    FWHM = []
    file = "Det_Theta.A"

    angle, intensity, error, count, _ = read_data_files(folder, file, D=1)
    
    # Plot Data
    fig, ax = plt.subplots()
    ax.set_xlabel("2$\Theta$ [deg]")
    ax.set_ylabel("Intensity")
    ax.errorbar(angle, intensity, yerr=error, fmt='o', color='b')
    plt.savefig(f"{folder}/distribution.png")
    plt.close()
    
    # Angles at which the intensity is measured
    angles = (np.linspace(10,90,17))
    angles_2 = (np.linspace(100,130,4))
    angles = np.concatenate((angles, angles_2))
    
    # Maximum width of the peak
    width = (np.ones(17)*5)
    width_2 = (np.ones(4)*10)
    widths = np.concatenate((width, width_2))
    
    for i in range(len(angles)):
        bin_min = angles[i] - widths[i]/2
        bin_max = angles[i] + widths[i]/2
        idx = np.where((angle > bin_min) & (angle < bin_max))

        intens_angles = intensity[idx]
        angles_angles = angle[idx]
        intensities.append(np.sum(intens_angles))
        FWHM.append(resolution(intens_angles, angles_angles)*2.355)
        
    
    intensities = np.array(intensities)
    if np.sum(intensities) == 0:
        er = 1
    FWHM = np.array(FWHM)
    return angles, intensities, FWHM, er

def get_offset(i, run, var, L1,L2,L3,mos1,mos2,col_h,col_v):
    folder = f"{main_folder}/{run}-{var}/{i}-offset"
    
    string = f"mcrun SLEIPNIR.instr -c --mpi=12 -n1e5 -d {folder} L0=2.6 L1={L1} L2={L2} L3={L3} mos1={mos1} mos2={mos2} col_h={col_h} col_v={col_v} d_spacing=3.355"
    
    call(string, shell=True)
    
    intensity, error, count, parameters = read_data_files(folder, "PSD_Final.dat", D=2)
    
    try:
        _, y_pos = average_pos(intensity, parameters)
        y_pos /= 100 # Convert to m
        er = 0
    except ZeroDivisionError:
        y_pos = 0
        er = 1
    
    return y_pos, er
    
def run_simulation(L0, d, d_spacing,n_var):
    var = d["Variable"]
    run = d["run"]
    
    os.makedirs(f"{main_folder}/{run}-{var}")
    
    vars = np.linspace(float(d[var][0]),float(d[var][1]),n_var,endpoint=True)
     
    for i in range(n_var):
        d[var] = vars[i]
        L1 = d["L1"]
        L2 = d["L2"]
        L3 = d["L3"]
        mos1 = d["mos1"]
        mos2 = d["mos2"]
        col_h = d["col_h"]
        col_v = d["col_v"]
        #HD = d["HD"]
        
        sample_y, er = get_offset(i, run, var, L1,L2,L3,mos1,mos2,col_h,col_v)
        
        worker(run,var,i,L0,L1,L2,L3,mos1,mos2,col_h,col_v, sample_y,d_spacing,er)
        

    d[var] = [vars[0],vars[-1]]
    write_path(d)
    administrator(run,var,n_var)
    return True

def worker(run,var,i,L0,L1,L2,L3,mos1,mos2,col_h, col_v, sample_y,d_spacing,er):
    
    # Add simulation to template
    if er==0:
        string = f"mcrun SLEIPNIR.instr -c --mpi=12 -n1e7 -d {main_folder}/{run}-{var}/{i} L0={L0} L1={L1} L2={L2} L3={L3} mos1={mos1} mos2={mos2} col_h={col_h} col_v={col_v} sample_y={sample_y} d_spacing={d_spacing}"
    else:
        string = f"mcrun SLEIPNIR.instr -c --mpi=12 -n1e5 -d {main_folder}/{run}-{var}/{i} L0={L0} L1={L1} L2={L2} L3={L3} mos1={mos1} mos2={mos2} col_h={col_h} col_v={col_v} d_spacing={d_spacing}"
    with open("template.scpt", "r") as f:
        lines = f.readlines()
    lines[5] = lines[5][:-1]
    lines[5] += f"worker_{order}_{run}_{i} \n"
    lines[21] = lines[21][:-1]
    lines[21] += f"{main_folder}/workers/worker_{run}_{i}.out \n"
    lines[22] = lines[22][:-1]
    lines[22] += f"{main_folder}/workers/worker_{run}_{i}.err \n"
    lines.append(string)
    
    with open(f"worker_{order}.scpt", "w") as f:    
        for line in lines:
            f.write(line)
    
    # Worker is called
    call(f"bsub < worker_{order}.scpt", shell=True)
    
    # Copy work script to script folder
    shutil.copy(f"worker_{order}.scpt", f"{main_folder}/scripts/{run}-{i}.scpt")
    
    # Wait 30 seconds
    time.sleep(30)
    return True

def administrator(run,var,n_var):
    # Check if all workers are done
    folder = f"{main_folder}/{run}-{var}"
    truth = np.zeros(n_var)
    j = 0
    while True:
        for i in range(n_var):
            if os.path.exists(f"{folder}/{i}/PSD_Final.dat"):
                truth[i] = 1
        if np.sum(truth) == n_var:
            break
        time.sleep(900)
        if j > 10:
            print("Truth: ", truth)
            raise ValueError("Workers are taking too long")
        j += 1
        if j % 3 == 0:
            # Check if any workers have failed
            for i in range(n_var):
                if truth[i] == 0:
                    print(f"Worker {i} failed")
                    try:
                        os.remove(f"{folder}/{i}")
                    except FileNotFoundError:
                        pass
                    shutil.copy(f"{main_folder}/scripts/{run}-{i}.scpt", f"worker_{order}.scpt")
                    call(f"bsub < worker_{order}.scpt", shell=True)
                    time.sleep(30)
    
    return 
            
        
def find_optimal(L0, d, d_spacing,n_var, prev_var):
    folder = f"{main_folder}/{d['run']}-{d['Variable']}"
    print(folder)
    vars = np.linspace(float(d[d['Variable']][0]),float(d[d['Variable']][1]),n_var,endpoint=True)
    I = np.zeros(n_var)
    R_RMS = np.zeros(n_var)
    for i in range(n_var):
        angles, intens, FWHM, error = read_data(f"{folder}/{i}", d)
        
        
        # Plot
        fig, ax = plt.subplots()
        ax.plot(angles, intens, color='b')
        ax.set_xlabel("2$\Theta$ [deg]")
        ax.set_ylabel("Intensity", color='b')
        
        ax2 = ax.twinx()
        ax2.plot(angles, FWHM, color='r')
        ax2.set_ylabel("FWHM", color='r')
        plt.savefig(f"{folder}/{i}/plot.png")
        plt.close()
        
        # Save data
        np.savetxt(f"{folder}/{i}/data.txt", np.array([angles, intens, FWHM]).T)
        if error == 1:
            I[i] = np.nan
            R_RMS[i] = np.nan
        intens = np.sum(intens)
        R_p = FWHM*L0**2*np.sin(angles*np.pi/180)
        R_RMS_i = np.sqrt(np.sum(R_p**2))
        I[i] = intens
        R_RMS[i] = R_RMS_i
        
    I = np.array(I)
    R_RMS = np.array(R_RMS)

    # Plot FOM
    for i in range(1,5):
        fig, ax = plt.subplots()
        FOM_temp = I/(R_RMS**i)
        popt, _ = curve_fit(third_order_poly, vars, FOM_temp)
        x = np.linspace(vars[0], vars[-1], 1000)
        FOM_fit = third_order_poly(x, *popt)
        ax.vlines(prev_var, min(FOM_temp), max(FOM_temp), color='k', linestyle='--')
        ax.plot(vars, FOM_temp, label="Data")
        ax.plot(x, FOM_fit, label="Fit")
        ax.set_xlabel("Variable")
        ax.set_ylabel("FOM")
        plt.savefig(f"{folder}/FOM_{i}.png")
        plt.close()

    # Find optimal
    FOM = I/(R_RMS**order)
    
    np.savetxt(f"{folder}/data.txt", np.array([vars, I, R_RMS ,FOM]).T)
    
    if prev_var < vars[0]:
        prev_var = vars[0]
    elif prev_var > vars[-1]:
        prev_var = vars[-1]
    
    # Gradient descent
    
    # Both FOM and vars are normalized
    norm = np.linspace(0,1,n_var)
    prev_norm = (prev_var-vars[0])/(vars[-1]-vars[0])
    min_max_FOM = (FOM-np.min(FOM))/(np.max(FOM)-np.min(FOM))
    
    # Find gradient
    popt, _ = curve_fit(third_order_poly, norm, min_max_FOM)
    grad = 3*popt[0]*prev_norm**2 + 2*popt[1]*prev_norm + popt[2]
    
    # Find optimal value
    opt = prev_norm + eta*grad    # Gradient descent, 1st order, for the maximum
    
    # Find closest value
    idx = np.argmin(np.abs(norm-opt))
    
    # Get new optimal value
    opt = vars[idx]
    FOM_val = FOM[idx]
    
    return opt, FOM_val


def third_order_poly(x, a,b,c,d):
    return a*x**3 + b*x**2 + c*x + d

def next_run(L0,d_spacing, n_var, FOM_ar):
    
    # Read previous
    d = parameters()
    var = d["Variable"]
    
    # Read previous value of variable
    prev_var = prev_parameters(var)
    
    # Find optimal value
    opt, FOM_val = find_optimal(L0, d, d_spacing,n_var, prev_var)
    FOM_ar.append(FOM_val)
    d[var] = opt
    
    print("Optimal value: ", d[var])

    # Pick next variable
    par_names = ["L1", "L2", "L3", "mos1", "mos2", "col_h","col_v"]
    idx = par_names.index(var)
    
    next_idx = idx + 1
    next_idx %= len(par_names)
 
    d["Variable"] = par_names[next_idx]
    
    # Check how many times this variable has been run
    with open(f'path_{order}.csv', 'r') as f:
        text = f.read()
    count = text.count(d["Variable"]) - 1
    
    d['run'] = count
    
    # Pick variable range
    if d["Variable"] == "L1":
        L2 = d["L2"]
        L1_max = (3.5/np.sin(Theta2(L0,d_spacing)))-L2
        L1_min = 0.5
        L1_mid = (L1_max-L1_min)/2
        L1_min_scaled = L1_mid + (L1_min-L1_mid)
        L1_max_scaled = L1_mid + (L1_max-L1_mid)
        d["L1"] = [L1_min_scaled, L1_max_scaled]
    elif d["Variable"] == "L2":
        L1 = d["L1"]
        L2_max = (3.5/np.sin(Theta2(L0,d_spacing)))-L1
        L2_min = 0.8
        L2_mid = (L2_max-L2_min)/2
        L2_min_scaled = L2_mid + (L2_min-L2_mid)
        L2_max_scaled = L2_mid + (L2_max-L2_mid)
        d["L2"] = [L2_min_scaled, L2_max_scaled]
    elif d["Variable"] == "L3":
        # d["L3"] = [0.75,2.25]
        L3_max = 2.25
        L3_min = 0.75
        L3_mid = (L3_max-L3_min)/2
        L3_min_scaled = L3_mid + (L3_min-L3_mid)
        L3_max_scaled = L3_mid + (L3_max-L3_mid)
        d["L3"] = [L3_min_scaled, L3_max_scaled]
    elif d["Variable"] == "mos1":
        # d["mos1"] = [30,90]
        mos1_max = 90
        mos1_min = 30
        mos1_mid = (mos1_max-mos1_min)/2
        mos1_min_scaled = mos1_mid + (mos1_min-mos1_mid)
        mos1_max_scaled = mos1_mid + (mos1_max-mos1_mid)
        d["mos1"] = [mos1_min_scaled, mos1_max_scaled]
    elif d["Variable"] == "mos2":
        # d["mos2"] = [30,90]
        mos2_max = 90
        mos2_min = 30
        mos2_mid = (mos2_max-mos2_min)/2
        mos2_min_scaled = mos2_mid + (mos2_min-mos2_mid)
        mos2_max_scaled = mos2_mid + (mos2_max-mos2_mid)
        d["mos2"] = [mos2_min_scaled, mos2_max_scaled]
        
    elif d["Variable"] == "col_h":
        # d["col_h"] = [10,60]
        col_h_max = 90
        col_h_min = 10
        col_h_mid = (col_h_max-col_h_min)/2
        col_h_min_scaled = col_h_mid + (col_h_min-col_h_mid)
        col_h_max_scaled = col_h_mid + (col_h_max-col_h_mid)
        d["col_h"] = [col_h_min_scaled, col_h_max_scaled]
        
        
    elif d["Variable"] == "col_v":
        # d["col_v"] = [10,60]
        col_v_max = 90
        col_v_min = 10
        col_v_mid = (col_v_max-col_v_min)/2
        col_v_min_scaled = col_v_mid + (col_v_min-col_v_mid)
        col_v_max_scaled = col_v_mid + (col_v_max-col_v_mid)
        d["col_v"] = [col_v_min_scaled, col_v_max_scaled]
    
    #elif d["Variable"] == "RV2":
    #    ttheta = Theta2(L0,d_spacing)
    #    RV_C = 2*d["L2"]/np.sin(ttheta)
    #    # d["RV2"] = [RV_C*0.75,RV_C*1.25]
    #    RV2_max = RV_C*1.25
    #    RV2_min = RV_C*0.75
    #    RV2_mid = (RV2_max-RV2_min)/2
    #    RV2_min_scaled = RV2_mid + (RV2_min-RV2_mid)
    #    RV2_max_scaled = RV2_mid + (RV2_max-RV2_mid)
    #    d["RV2"] = [RV2_min_scaled, RV2_max_scaled]
    
    #elif d["Variable"] == "HD":
    #    d["HD"] = [0.1,0.5]
    #    HD_max = 0.5
    #    HD_min = 0.1
    #    HD_mid = (HD_max-HD_min)/2
    #    HD_min_scaled = HD_mid + (HD_min-HD_mid)
    #    HD_max_scaled = HD_mid + (HD_max-HD_mid)
    #    d["HD"] = [HD_min_scaled, HD_max_scaled]
        
    
    return d, FOM_ar

def cycle(L0,d_spacing,n_var, n_cycles=0):
    FOM_ar = []
    os.makedirs(f"{main_folder}/plots", exist_ok=True)
    if os.path.exists(f"{main_folder}/plots/data.txt"):
        with open(f"{main_folder}/plots/data.txt", "r") as f:
            lines = f.readlines()
        for line in lines:
            FOM_ar.append(float(line))        
    d, FOM_ar = next_run(L0, d_spacing, n_var, FOM_ar)
    with open(f"{main_folder}/plots/data.txt", "w") as f:
        for FOM in FOM_ar:
            f.write(str(FOM)+"\n")
    t0 = time.time()
    tmax = 60*60*24*3
    t_iter_array = []
    t_iter=0
    
    
    with open(f'path_{order}.csv', 'r') as f:
        lines = f.readlines()
        i = len(lines)
    
    if n_cycles == 0:
        while time.time() - t0 < tmax-2*t_iter:
            print("-"*len(f"Cycle {i}"))
            print(f"Cycle {i} - Variable: {d['Variable']}")
            print("-"*len(f"Cycle {i}"))
            t_start = time.time()
            run_simulation(L0, d, d_spacing,n_var)
            d, FOM_ar = next_run(L0, d_spacing, n_var, FOM_ar)
            
            # Plot FOM pr. cycle
            fig, ax = plt.subplots()
            ax.plot(FOM_ar)
            ax.set_xlabel("Cycle")
            ax.set_ylabel("FOM")
            plt.savefig(f"{main_folder}/plots/FOM_{i}.png")
            plt.close()
            
            with open(f"{main_folder}/plots/data.txt", "w") as f:
                for FOM in FOM_ar:
                    f.write(str(FOM)+"\n")
            
            t_end = time.time()
            print("-"*len(f"Cycle done"))
            print("Cycle done - Variable: ", d["Variable"])        
            print("-"*len(f"Cycle done"))    
            t_iter_array.append(t_end-t_start)
            t_iter = np.mean(t_iter_array)
            i += 1
    else:
        for j in range(n_cycles):
            print("-"*len(f"Cycle {i}"))
            print(f"Cycle {i} - Variable: {d['Variable']}")
            print("-"*len(f"Cycle {i}"))
            t_start = time.time()
            run_simulation(L0, d, d_spacing,n_var)
            d = next_run(L0, d_spacing, n_var)
            t_end = time.time()
            print("-"*len(f"Cycle done"))
            print("Cycle done - Variable: ", d["Variable"])        
            print("-"*len(f"Cycle done"))    
            t_iter_array.append(t_end-t_start)
            t_iter = np.mean(t_iter_array)
            i += 1
        i = n_cycles
    print("Time: ", time.time()-t0)
    print("Time per iteration: ", t_iter)
    print("Iterations: ", i)

def first_cycle(L0,d_spacing,n_var,d1):
    os.makedirs(main_folder, exist_ok=True)
    os.makedirs(main_folder+"/workers", exist_ok=True)
    os.makedirs(main_folder+"/scripts", exist_ok=True)
    print("-"*len(f"Cycle 0"))
    print(f"Cycle 0 - Variable: {d1['Variable']}")
    print("-"*len(f"Cycle 0"))
    run_simulation(L0, d1, d_spacing,n_var)
    print("-"*len(f"Cycle done"))
    print("Cycle done - Variable: ", d1["Variable"])
    print("-"*len(f"Cycle done"))
    return True

d_SI111 = 3.355
L0 = 2.6
order = 4
main_folder = f"L0-{L0}-{order}"
eta = 0.1


n_var = 41
# Cycle 1
#d = {"run": 0, "L1": [1,(3.5/np.sin(Theta2(L0,d_SI111)))-1.5], "L2": 1.5, "L3": 1.5, "mos1": 60, "mos2": 60, "col_h": 60, "col_v": 60, "Variable": "L1"}
#first_cycle(L0,d_SI111,n_var,d)
# Cycles Next
cycle(L0, d_SI111,n_var=n_var)

