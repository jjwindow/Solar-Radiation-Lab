"""
Jack J. Window
3rd Year Physics Undergraduate
Blackett Laboratory, Imperial College London

07.04.2020

Data analysis script for Solar Radiation 3rd year experiment.
Data taken on 09.03.2020 during cloudy conditions to study the probability distribution
of cloud optical depths in London.
"""

# ----------------------------------- IMPORT MODULES -----------------------------------
import csv
from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt
from math import floor
from scipy.special import gamma, wofz
from scipy.stats import beta
from plotSetup import plotSetup

# --------------------------------------------------------------------------------------

def readData(filename):
    """
    Function to read data from a pre-organised CSV of solar instrument measurements.
    
    CSV is assumed to be of the format ->

    NUM | LIGHTMETER [KLux] | SOLARIMETER [mV]  <- Headers here given as indication. 
    ----------------------------------             CSV must have NO HEADERS to process
     1  |        lm1        |      sol1             the data properly.
     .  |        .          |       .
     .  |        .          |       .
     .  |        .          |       .
     60 |        lmN        |       solN           <- For N data points
    ----------------------------------

    PARAMETERS
    -------------------------------------------------------------------------------------
    filename    :   str, the name of the csv from which to read data if stored in the 
                    same directory as the script. Otherwise this should be a full filepath.
    
    RETURNS
    -------------------------------------------------------------------------------------
    N           :   int, number of measurements in csv.
    raw_lm      :   list, contains all raw lightmeter measurements (no conversion factor)
    raw_sol     :   list, contails all raw solarimeter measurements (no conversion factor)
    """   

    # read raw measurements from CSV   
    raw_lm = []
    raw_sol = []
    with open(filename) as solar:
        rawData = csv.reader(solar, delimiter=',')
        for row in rawData:
            # Populate column arrays with data from CSV, using format in docstring
            raw_lm.append(row[1])
            raw_sol.append(row[2])

    N = len(raw_lm)

    # Perform conversions
    # Solarimeter: mV -> W/m^2 = x14.98
    # Lightmeter: Klux -> flux = x7.9 (0.0079 flux per lux)
    lm = [7.9*float(_lm)/1000 for _lm in raw_lm]
    sol = [14.98*float(_sol) for _sol in raw_sol]

    return N, lm, sol

def clearSkyFlux(t):
    """
    take t as decimal, return clear sky model flux in W/m^2.
    """
    model_flux = {10.9 : 454, 11.0 : 465, 11.1 : 476, 11.2 : 486, 11.3 : 495, 11.4 : 504, 11.5 : 513, 11.6 : 521, 11.7 : 529, 11.8 : 536, 11.9 : 543, 12.0 : 549}

    if (t < 10.9) or (t>11.9):
        raise Exception("Time not in measured values.")

    if t in model_flux:
        return model_flux[t]
    else:
        # Interpolate
        t_below = round(floor(100*t)/100, 1)
        t_above = round(t_below + 0.1, 1)
        f_below = model_flux[t_below]
        f_above = model_flux[t_above]
        m = (f_above-f_below)/(t_above-t_below)
        return f_below + (t-t_below)*m

# # Clear sky solar flux over time of data taking.
# t_plt = np.arange(10.9, 11.9, 0.05)
# f_plt = [clearSkyFlux(t) for t in t_plt]
# plt.plot(t_plt, f_plt)
# plt.show()

def transmittance(t, x):
    """
    Returns the transmittance from a value of measured solar flux in a cloudy sky
    """
    # Clear sky solar flux
    f_cs = clearSkyFlux(t)
    T = x/f_cs
    # Measurement error est. at 2% due to tilt
    sig_T = 0.02*T
    return T, sig_T

def opticalDepth(T, sig_T):
    albedo = 0.265 # Charlotte's value for asphalt using solarimeter.
    sig_a = 0.027
    # formula from Si-Woo's report
    u = (1/T - 1)*(1/(0.1153-0.1226*albedo))
    # constants in eq.
    c = 0.1153
    d = 0.1126
    sig_u = np.sqrt((sig_T/(c-d*albedo))**2 / T**4 + (sig_a*d*(1/T - 1))**2 / (c-d*albedo)**4)
    return u, sig_u

def opticalDepthPlot(plot = True):
    filename = 'cloud_distribution.csv'
    N, lm, sol = readData(filename)
    time = [10.9 + i/60 for i in range(N)]
    # Array of measurement-error pairs for transmittance
    T_with_sig = [transmittance(t, f_cloud) for t, f_cloud in zip(time, sol)]
    # Transposed list
    # Tvals_all_trans = [list(T) for T in zip(*Tvals_all)]
    # T_arr = Tvals_all_trans[0]  # Transmittance values
    # T_sig = Tvals_all_trans[1]  # Absolute uncertainties
    u_arr = [opticalDepth(T, sig) for T, sig in T_with_sig]
    u_arr_trans = [list(u) for u in zip(*u_arr)]
    u = u_arr_trans[0]
    sig_u = u_arr_trans[1]
    if plot:
        font_prop = plotSetup('Time of Day (hours)', 'Optical Depth', xscale='linear', yscale = 'linear')
        plt.errorbar(time, u, yerr=sig_u, fmt='rx')
        plt.plot(time, u, '-k')
        plt.show()
    return time, u, sig_u

# opticalDepthPlot()

def normal(x, m, s):
    return (1/(np.sqrt(2*np.pi)*s))*np.exp(-0.5*((x-m)**2)/s**2)

def poisson(x, l):
    if type(x) is list or type(x) is np.ndarray:
        return [(l**_x)*np.exp(-l)/gamma(_x+1) for _x in x]
    else:
        return (l**x)*np.exp(-l)/gamma(x+1)

# def G(x, alpha):
#     """ Return Gaussian line shape at x with HWHM alpha """
#     return np.sqrt(np.log(2) / np.pi) / alpha\
#                              * np.exp(-(x / alpha)**2 * np.log(2))

def L(x, m, gamma):
    """ Return Lorentzian line shape at x with HWHM gamma """
    return gamma / np.pi / ((x-m)**2 + gamma**2)

def V(x, m, sigma, gamma):
    """
    Return the Voigt line shape at x with Lorentzian component HWHM gamma
    and Gaussian component HWHM alpha.

    """
    # sigma = alpha / np.sqrt(2 * np.log(2))

    return np.real(wofz(((x-m) + 1j*gamma)/sigma/np.sqrt(2))) / sigma\
                                                           /np.sqrt(2*np.pi)

def GplusL(x ,m_g, m_l, sigma, gamma):
    return normal(x, m_g, sigma) + L(x, m_l, gamma)

def plotProbDist(fitnorm = False, fitpoisson = False, fitall = False, fitgl = True, fitgplusl = True):
    time, u, sig_u = opticalDepthPlot(False)
    font_prop = plotSetup('Optical Depth', 'Probability density', xscale='linear', yscale='linear')
    # plt.hist(u, bins = 'auto', density=True, histtype='stepfilled', color = 'k')
    p, edges = np.histogram(u, bins='auto', density=True)
    u_bins = [0.5*(a+b) for a, b in zip(edges, edges[1:])]
    plt.plot(u_bins, p, 'rx', label='Measured Data', markersize = 10)
    fit_vals = {}

    def fitG(bins, p):
        p0 = [80, 10]
        n_popt, n_pcov = curve_fit(normal, bins, p, p0=p0)
        n_perr = np.sqrt(np.diag(n_pcov))
        x = np.linspace(min(bins)-10, max(bins)+20, 1000)
        plt.plot(x, [normal(_x, n_popt[0], n_popt[1]) for _x in x], '--k', label = "Gaussian")
        fit_vals["norm"] = list(zip(n_popt, n_perr))
        return list(zip(n_popt, n_perr))

    def fitP(bins, p):
        p0 = [60]
        p_popt, p_pcov = curve_fit(poisson, bins, p, p0=p0)
        p_perr = np.sqrt(np.diag(p_pcov))
        x = np.linspace(min(bins)-10, max(bins)+20, 1000)
        plt.plot(x, [poisson(_x, p_popt[0]) for _x in x], '--k', label = "Poisson")
        fit_vals["Poisson"] = list(zip(p_popt, p_perr))
        return list(zip(p_popt, p_perr))

    def fitV(bins, p):
        p0 = [60, 20, 10]
        v_popt, v_pcov = curve_fit(V, bins, p, p0=p0)
        v_perr = np.sqrt(np.diag(v_pcov))
        x = np.linspace(min(bins)-10, max(bins)+20, 1000)
        plt.plot(x, V(x, v_popt[0], v_popt[1], v_popt[2]), '--k', label = "Voigt")
        fit_vals["Voigt"] = list(zip(v_popt, v_perr))
        return list(zip(v_popt, v_perr))

    def fitL(bins, p):
        p0 = [60, 5]
        l_popt, l_pcov = curve_fit(L, bins, p, p0=p0)
        l_perr = np.sqrt(np.diag(l_pcov))
        x = np.linspace(min(bins)-10, max(bins)+20, 1000)
        plt.plot(x, L(x, l_popt[0], l_popt[1]), ':', label = "Lorentz")
        fit_vals["Lorentz"] = list(zip(l_popt, l_perr))
        return list(zip(l_popt, l_perr))

    def fitGplusL():
        p0 = [80, 60, 20, 5]
        gl_popt, gl_pcov= curve_fit(GplusL, u_bins, p, p0=p0)
        gl_perr = np.sqrt(np.diag(gl_pcov))
        x = np.linspace(min(u_bins)-10, max(u_bins)+20, 1000)
        plt.plot(x, GplusL(x, gl_popt[0], gl_popt[1], gl_popt[2], gl_popt[3]), 'k-', label = "Gauss + Lorentz")
        fit_vals["G+L"] = list(zip(gl_popt, gl_perr))
        return list(zip(gl_popt, gl_perr))
    
    def fitGL():
        l_slice = [0, 1, 2, 6, 7]
        g_slice = [0, 2, 3, 4, 5, 6, 7]
        g_bins = [u_bins[g] for g in g_slice]
        p_g = [p[g] for g in g_slice]
        l_bins = [u_bins[l] for l in l_slice]
        p_l = [p[l] for l in l_slice]
        fitL(l_bins, p_l)
        fitG(g_bins, p_g)

    if fitgl:
        params = fitGL()
        
    if fitgplusl:
        params = fitGplusL()
        print(f"Gaussian Parameters\nMean +/- error: {params[0]}\nStdev +/- error: {params[2]}")
        print(f"Lorenzian Parameters\nMean +/- error: {params[1]}\nGamma +/- error: {params[3]}")
        x = np.linspace(min(u_bins)-10, max(u_bins)+20, 1000)
        plt.plot(x, [normal(_x, params[0][0], params[2][0]) for _x in x], '--k', label = "Gaussian")
        plt.plot(x, L(x, params[1][0], params[3][0]), ':k', label = "Lorentz")

    if fitall:
        fitG(u_bins, p)
        fitV(u_bins, p)
        fitL(u_bins, p)
    else:
        if fitnorm:
            fitG(u_bins, p)
        if fitpoisson:
            fitP(u_bins, p)

    plt.legend(prop=font_prop)
    plt.show()
    return fit_vals

plotProbDist(fitgl=False)

