"""
Functions for use in McStas data analysis. Import this module to use the functions in your script.

Written by: Nicolai Amin

Date: 2024-06-06

Functions:
- read_data_files(folder_path, file_name, D=2)
- xy(parameters, intensity)
- plot_heatmap(data, parameters, logarithm=False, with_fit=False)
- gaussian_2d(xy, amplitude, xo, yo, sigma_x, sigma_y, theta)
- fit_gaussian(data, parameters)
- pearsons_chi_square_test(data, fit_data)
- r_squared(data, fit_data)
- mccode_int(folder, monitor="l_monitor_0_I")
- wstd(values, weights)
- average_width(intensity, parameters)
- average_pos(intensity, parameters)
- resolution(intensity, x)
- E_resolution(intensity, x)
- FWHM_calc(x, intensity)
- bragg_theta(f, d)
- FWHM_angle_fit(theta, U, V, W)
- FWHM_fitting(angle, FWHM)
- D1B_plot(ax)
- ThetaTOF(folder)
"""

import os
import math
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
# Set the matplotlib parameters for LaTeX rendering
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.sans-serif": ["Computer Modern Roman"],
    "font.size": 14,
    "font.weight": "bold",
    "legend.fontsize": 12,
    "legend.edgecolor": [0.2, 0.2, 0.2],
    "axes.linewidth": 1.75,
    "axes.titlesize": 20,
    'text.latex.preamble': r'\boldmath',
    "figure.autolayout": True})
from scipy.optimize import curve_fit
import scipy.stats as stats
import types

def bottom_offset(self, bboxes, bboxes2):
    """
    Sets the position of the offset text at the bottom of the plot.

    Parameters:
    - bboxes: The bounding boxes of the plot elements.
    - bboxes2: The second set of bounding boxes.

    Returns:
    None
    """
    bottom = self.axes.bbox.ymax
    self.offsetText.set(va="top", ha="left")
    self.offsetText.set_position(
            (-0.1, bottom + 5*self.OFFSETTEXTPAD * self.figure.dpi / 72.0))
def register_bottom_offset(axis, func):
    """
    Register a custom function to update the offset text position for the given axis.

    Parameters:
    - axis: The axis object for which the offset text position needs to be updated.
    - func: The custom function that will be used to update the offset text position.

    Returns:
    None
    """
    axis._update_offset_text_position = types.MethodType(func, axis)

def bragg(d, L):
    """
    Calculate the Bragg angle for a given lattice spacing and wavelength.

    Parameters:
    d (float): Lattice spacing.
    L (float): Wavelength.

    Returns:
    float: Bragg angle in radians.
    """
    return 2 * np.arcsin(L / (2 * d))

# Read the data files
def read_data_files(folder_path, file_name, D=2):
    """
    Read data files from a specified folder path and extract the data based on the given file name and dimension.

    Parameters:
        folder_path (str): The path to the folder containing the data files.
        file_name (str): The name of the file to be read.
        D (int, optional): The dimension of the data. Default is 2.

    Returns:
        tuple: A tuple containing the extracted data and parameters. The structure of the tuple depends on the dimension (D) parameter.
            - If D=2, the tuple contains:
                - intensity (numpy.ndarray): The intensity data.
                - error (numpy.ndarray): The error data.
                - count (numpy.ndarray): The count data.
                - parameters (dict): A dictionary containing the extracted parameters.
            - If D=1, the tuple contains:
                - lam (numpy.ndarray): The wavelength data.
                - intensity (numpy.ndarray): The intensity data.
                - error (numpy.ndarray): The error data.
                - count (numpy.ndarray): The count data.
                - parameters (dict): A dictionary containing the extracted parameters.
    """
    parameters = {}
    for root, dirs, files in os.walk(folder_path):
        for file in files:
            if file == file_name and D==2:
                file_path = os.path.join(root, file)
                with open(file_path, 'r') as f:
                    lines = f.readlines()
                    lines = [line.strip() for line in lines]
                    for line in lines:
                        if line.startswith("#"):
                            if line.startswith("# Param"):
                                key, value = line.split("=")[0].strip().split(" ")[-1], line.split("=")[1].strip()
                                parameters[key] = value
                            else:
                                key, value = line.split(":")[0].strip().split(" ")[-1], line.split(":")[1].strip()
                                parameters[key] = value
                    lines = [list(map(float, line.split())) for line in lines if not line.startswith("#")]
                    data = np.array(lines)
                    
                    # Separate the data into three chunks
                    intensity = data[0:data.shape[0]//3, :]
                    error = data[data.shape[0]//3:2*(data.shape[0]//3), :]
                    count = data[2*(data.shape[0]//3):, :]
                    return intensity, error, count, parameters
            if file == file_name and D == 1:
                file_path = os.path.join(root, file)
                with open(file_path, 'r') as f:
                    lines = f.readlines()
                    lines = [line.strip() for line in lines]
                    for line in lines:
                        if line.startswith("#"):
                            if line.startswith("# Param"):
                                key, value = line.split("=")[0].strip().split(" ")[-1], line.split("=")[1].strip()
                                parameters[key] = value
                            else:
                                key, value = line.split(":")[0].strip().split(" ")[-1], line.split(":")[1].strip()
                                parameters[key] = value
                    lines = [list(map(float, line.split())) for line in lines if not line.startswith("#")]
                    data = np.array(lines)
                    lam = data[:, 0]
                    intensity = data[:, 1]
                    error = data[:, 2]
                    count = data[:, 3]
                    return lam, intensity, error, count, parameters
                    
    

# Define the x and y meshgrid
def xy(parameters, intensity):
    """
    Generate x and y coordinates based on the given parameters and intensity shape.

    Parameters:
        parameters (dict): A dictionary containing the xylimits as a string in the format "x1 x2 y1 y2".
        intensity (ndarray): An array representing the intensity.

    Returns:
        x (ndarray): An array of x-coordinates.
        y (ndarray): An array of y-coordinates.
    """
    xylimits = parameters["xylimits"].split(" ")
    x = np.linspace(float(xylimits[0]), float(xylimits[1]), intensity.shape[1])
    y = np.linspace(float(xylimits[2]), float(xylimits[3]), intensity.shape[0])
    x, y = np.meshgrid(x, y)
    return x, y

# Plot the heatmap
def plot_heatmap(data, parameters, logarithm=False, with_fit=False):
    """
    Plot a heatmap of the given data.

    Parameters:
    - data: numpy.ndarray
        The 2D array of data to be plotted.
    - parameters: dict
        A dictionary containing the plot parameters.
        It should have the following keys:
        - 'xlabel': str
            The label for the x-axis.
        - 'ylabel': str
            The label for the y-axis.
        - 'zlabel': str
            The label for the colorbar.
        - 'title': str
            The title of the plot.
    - logarithm: bool, optional
        If True, the data will be plotted in logarithmic scale.
        Default is False.
    - with_fit: bool, optional
        If True, the plot will include a fitted Gaussian curve.
        Default is False.
    """
    x, y = xy(parameters,data)
    if logarithm:
        data = np.log(data)
    # Plot the heatmap with fit
    if with_fit:
        plt.figure()
        plt.subplot(1, 2, 1)
        plt.imshow(data, cmap='hot', extent=[x.min(), x.max(), y.min(), y.max()])
        plt.xlabel(parameters['xlabel'])
        plt.ylabel(parameters['ylabel'])
        plt.title("Data")
        plt.colorbar(label=parameters['zlabel'])  # Add the colorbar label
        plt.subplot(1, 2, 2)
        params = fit_gaussian(data,parameters)
        fit_data = gaussian_2d((x, y), *params).reshape(data.shape)
        plt.imshow(fit_data, cmap='hot', extent=[x.min(), x.max(), y.min(), y.max()])
        plt.xlabel(parameters['xlabel'])
        plt.ylabel(parameters['ylabel'])
        plt.title("Fit")
        plt.colorbar(label=parameters['zlabel'])  # Add the colorbar label
        plt.show()
    # Plot the heatmap without fit
    else:
        plt.imshow(data, cmap='hot', extent=[x.min(), x.max(), y.min(), y.max()])
        plt.xlabel(parameters['xlabel'])
        plt.ylabel(parameters['ylabel'])
        plt.title(parameters['title'])
        plt.colorbar(label=parameters['zlabel'])  # Add the colorbar label
        plt.show()

# Define the 2D Gaussian function
def gaussian_2d(xy, amplitude, xo, yo, sigma_x, sigma_y, theta):
    """
    Calculate the value of a 2D Gaussian function at a given point.

    Parameters:
    xy (tuple): The coordinates (x, y) at which to evaluate the Gaussian function.
    amplitude (float): The amplitude of the Gaussian function.
    xo (float): The x-coordinate of the center of the Gaussian function.
    yo (float): The y-coordinate of the center of the Gaussian function.
    sigma_x (float): The standard deviation of the Gaussian function in the x-direction.
    sigma_y (float): The standard deviation of the Gaussian function in the y-direction.
    theta (float): The rotation angle of the Gaussian function.

    Returns:
    float: The value of the Gaussian function at the given coordinates.
    """
    x, y = xy
    xo = float(xo)
    yo = float(yo)
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = amplitude * np.exp(- (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) + c*((y-yo)**2)))
    return g.ravel()

# Fit the 2D Gaussian function
def fit_gaussian(data, parameters):
    """
    Fits a 2D Gaussian function to the given data.

    Parameters:
    - data: 2D array-like object representing the data to fit the Gaussian to.
    - parameters: Additional parameters needed for the fitting process.

    Returns:
    - popt: Optimal values for the parameters of the Gaussian function.

    """
    x, y = xy(parameters, data)
    initial_guess = [np.max(data), np.mean(x), np.mean(y), 1, 1, 0]  # Initial guess for the parameters
    popt, pcov = curve_fit(gaussian_2d, (x, y), data.ravel(), p0=initial_guess)  # Perform the curve fitting
    return popt

# Calculate the Pearson's chi-square test
def pearsons_chi_square_test(data, fit_data):
    """
    Perform Pearson's chi-square test of independence.

    Parameters:
    - data (ndarray): Observed data.
    - fit_data (ndarray): Expected data.

    Returns:
    - chi2 (float): The chi-square test statistic.
    - p_value (float): The p-value associated with the test statistic.
    """
    observed = data.ravel()
    expected = fit_data.ravel()
    chi2, p_value = stats.chisquare(f_obs=observed, f_exp=expected)
    return chi2, p_value

# Calculate the R squared value
def r_squared(data, fit_data):
    """
    Calculate the coefficient of determination (R-squared) for a given data set and its corresponding fit data.

    Parameters:
    data (array-like): The actual data.
    fit_data (array-like): The fitted data.

    Returns:
    float: The coefficient of determination (R-squared) value.
    """
    y_mean = np.mean(data)
    ss_total = np.sum((data - y_mean)**2)
    ss_residual = np.sum((data - fit_data)**2)
    r2 = 1 - (ss_residual / ss_total)
    return r2

def mccode_int(folder, monitor="l_monitor_0_I"):
    """
    Read and extract data from a mccode.dat file.

    Parameters:
    folder (str): The path to the folder containing the mccode.dat file.
    monitor (str, optional): The name of the monitor to extract data from. Defaults to "l_monitor_0_I".

    Returns:
    tuple: A tuple containing two numpy arrays, x and y, representing the x and y coordinates of the monitor data.
    """
    file_path = os.path.join(folder, "mccode.dat")
    with open(file_path, 'r') as f:
        lines = f.readlines()
        lines = [line.strip() for line in lines] 
        for i in lines:
            if i.startswith("# variables"):
                monitors = i.split(" ")
                j = monitors.index(monitor)-2
        lines = [list(map(float, line.split())) for line in lines if not line.startswith("#")]
        data = np.array(lines)
    x = data[:, 0]
    y = data[:, j]
    return x, y

def wstd(values, weights):
    """
    Return the weighted average and standard deviation.

    They weights are in effect first normalized so that they 
    sum to 1 (and so they must not all be 0).

    values, weights -- NumPy ndarrays with the same shape.
    """
    try:
        average = np.average(values, weights=weights)
    except ZeroDivisionError:
        return np.nan
    # Fast and numerically precise:
    variance = np.average((values-average)**2, weights=weights)
    return math.sqrt(variance)

def average_width(intensity, parameters):
    """
    Calculate the average width of the intensity distribution along the x and y axes.

    Parameters:
    intensity (ndarray): 2D array representing the intensity distribution.
    parameters (dict): Dictionary containing the parameters for calculating the width.

    Returns:
    x_std (float): Standard deviation of the intensity distribution along the x axis.
    y_std (float): Standard deviation of the intensity distribution along the y axis.
    """
    xylimits = parameters["xylimits"].split(" ")
    x = np.linspace(float(xylimits[0]), float(xylimits[1]), intensity.shape[1])
    y = np.linspace(float(xylimits[2]), float(xylimits[3]), intensity.shape[0])
    x_hist = np.sum(intensity, axis=0)
    y_hist = np.sum(intensity, axis=1)
    x_std = wstd(x, x_hist)
    y_std = wstd(y, y_hist)
    return x_std, y_std

def average_pos(intensity, parameters):
    """
    Calculate the average position of intensity values in a 2D array.

    Parameters:
    intensity (ndarray): 2D array of intensity values.
    parameters (dict): Dictionary containing parameters for calculating the average position.
                       It should contain the key "xylimits" which is a string representing the limits of x and y values.

    Returns:
    tuple: A tuple containing the average x position and average y position.
    """
    xylimits = parameters["xylimits"].split(" ")
    x = np.linspace(float(xylimits[0]), float(xylimits[1]), intensity.shape[1])
    y = np.linspace(float(xylimits[2]), float(xylimits[3]), intensity.shape[0])
    x_hist = np.sum(intensity, axis=0)
    y_hist = np.sum(intensity, axis=1)
    x_pos = np.average(x, weights=x_hist)
    y_pos = np.average(y, weights=y_hist)
    return x_pos, y_pos

def resolution(intensity, x):
    """
    Calculate the resolution of the given intensity and x values.

    Parameters:
    intensity (float): The intensity value.
    x (float): The x value.

    Returns:
    float: The calculated resolution.
    """
    x_std = wstd(x, intensity)
    return x_std

def E_resolution(intensity, x):
    """
    Calculate the energy resolution.

    Parameters:
    intensity (float): The intensity value.
    x (float): The x value.

    Returns:
    float: The energy resolution.
    """
    E = 81.82 / (x ** 2)
    E_std = wstd(E, intensity)
    return E_std

def FWHM_calc(x, intensity):
    """
    Calculate the full width at half maximum (FWHM) of a peak.

    Parameters:
    x (array-like): The x-values of the peak.
    intensity (array-like): The intensity values of the peak.

    Returns:
    float: The FWHM of the peak.
    """
    half_max = np.max(intensity) / 2
    left_idx = np.where(intensity >= half_max)[0][0]
    right_idx = np.where(intensity >= half_max)[0][-1]
    FWHM = x[right_idx] - x[left_idx]
    return FWHM

def bragg_theta(f, d):
    """
    Calculate the Bragg angle (in radians) for a given frequency and lattice spacing.
    
    Parameters:
        f (float): The frequency of the wave.
        d (float): The lattice spacing.
    
    Returns:
        float: The Bragg angle in radians.
    """
    return 2 * np.arcsin(f / (2 * d))

def FWHM_angle_fit(theta, U, V, W):
    """
    Calculate the Full Width at Half Maximum (FWHM) angle fit.

    Parameters:
    theta (float): The angle in degrees.
    U (float): Coefficient U.
    V (float): Coefficient V.
    W (float): Coefficient W.

    Returns:
    float: The FWHM angle fit.

    """
    return np.sqrt(W + V*(np.tan(np.deg2rad(theta/2))) + U*(np.tan(np.deg2rad(theta/2)))**2)

def FWHM_fitting(angle, FWHM):
    """
    Fits the FWHM (Full Width at Half Maximum) values to the given angles using curve fitting.

    Parameters:
    angle (array-like): The array of angles.
    FWHM (array-like): The array of FWHM values.

    Returns:
    tuple: A tuple containing the optimized parameters for the curve fit.
    """
    popt, pcov = curve_fit(FWHM_angle_fit, angle, FWHM, maxfev=10000, p0=[0.01, -0.1, 0.5], bounds=([-np.inf, -np.inf, 0], [np.inf, 0, np.inf]))
    return popt


def D1B_plot(ax):
    """
    Plot the D1B data on the given axes.

    Parameters:
    - ax (matplotlib.axes.Axes): The axes on which to plot the data.

    Returns:
    None
    """
    
    folder_path = 'D1B'
    save_folder = "Plots"

    save_folder = os.path.join(folder_path, save_folder)
    with open("{}/data.txt".format(save_folder), "r") as file:
        lines = file.readlines()
    deg = []
    FWHM = []
    if len(lines) > 0:
        i = 0
        while i < len(lines):
            if lines[i].startswith("FWHM"):
                FWHM_temp = []
                deg_temp = []
                i += 1
                while i < len(lines) and lines[i] != "\n":
                    if "D1B" in folder_path:
                        deg_temp.append(-float(lines[i].split(",")[0]))
                    else:
                        deg_temp.append(float(lines[i].split(",")[0]))
                    FWHM_temp.append(float(lines[i].split(",")[1]))
                    i += 1
                deg.append(deg_temp)
                FWHM.append(FWHM_temp)
            i += 1
    deg = np.array(deg[0])*-1
    FWHM  = np.array(FWHM[0])   

    a = np.where(deg - 10 > 0)[0]
    #print(a)
    #print(deg[a], FWHM[a])
    ax.scatter(deg[a], FWHM[a], label="D1B", marker='x', color="black")
    popt = FWHM_fitting(deg[a],FWHM[a])
    x = np.linspace(0, 130, 1000)
    y = FWHM_angle_fit(x, *popt)
    ax.plot(x, y, label="Fit, D1B, U={:.2f}, V={:.2f}, W={:.2f}".format(*popt), color="black")
    
def ThetaTOF(folder):
    """
    Calculate the FWHM (Full Width at Half Maximum) and degree values for each time-of-flight (ToF) data point.
    
    Parameters:
    folder (str): The path to the folder containing the data files.
    
    Returns:
    deg (list): List of degree values for each ToF data point.
    FWHM_ToF (list): List of FWHM values for each ToF data point.
    t (list): List of time values for each ToF data point.
    """
    
    for i in os.listdir(folder):
        if i.endswith(".A_t"):
            file = i
            break
    intensity, error, count, parameters = read_data_files(folder, file, D=2)
    x, y = xy(parameters, intensity)
    FWHM_ToF = []
    deg = []
    t = []
    for i in range(len(intensity[:,0])):
        if np.sum(intensity[i,:]) < 10:
            continue
        else:
            t.append(y[i][0])
            FWHM = resolution(intensity[i,:],x[0])*2.355
            deg_t = np.average(x[0],weights=intensity[i,:])
            FWHM_ToF.append(FWHM)
            deg.append(deg_t)
    return deg, FWHM_ToF, t
    