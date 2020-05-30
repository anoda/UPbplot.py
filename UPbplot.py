#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# This is a script for calculation and visualization tool of U-Pb age
# data.  The script was written in Python 3.6.6

# Last updated: 2020/05/30 16:46:56.
# Written by Atsushi Noda
# License: Apache License, Version 2.0

# This software is provided "as is", without any warranty or guarantee
# of its usability or fitness for any purpose.  I don't provide
# support on the usage of the included software.  See the
# documentation and help resources on their websites if you need such
# help.

# __version__ = "0.0.1"           # Apr/01/2016
# __version__ = "0.0.2"           # Sep/23/2016
# __version__ = "0.0.3"           # Nov/15/2016
# __version__ = "0.0.4"             # Mar/01/2017
# __version__ = "0.0.5"             # Oct/30/2017
# __version__ = "0.0.6"  # Sep/12/2018
# __version__ = "0.0.7"  # Jan/24/2019
# __version__ = "0.0.8"  # Jun/06/2019
# __version__ = "0.0.9"  # Jul/03/2019
# __version__ = "0.1.0"  # May/22/2020
# __version__ = "0.1.1"  # May/23/2020
# __version__ = "0.1.2"  # May/24/2020
# __version__ = "0.1.3"  # May/27/2020
__version__ = "0.1.4"  # May/30/2020

# [Citation]
#
# Noda, Atsushi (2016) UPbplot.py: A python script for U-Pb age data
# analysis. Open-File Report, no. 634, Geological Survey of Japan,
# AIST.
#
# Noda, Atsushi (2017) A new tool for calculation and visualization of
# U--Pb age data: UPbplot.py.  Bulletin of the Geological Survey of
# Japan, vol. 68(3), p. 131-140, doi: 10.9795/bullgsj.68.131.

# [Preparation]
#
# 1. When you use this script (UPbplot.py) from the source code,
#    additional libraries of
#
#    Numpy: http://www.numpy.org
#    Matplotlib: http://matplotlib.org
#    pandas: http://pandas.pydata.org
#    SciPy: https://www.scipy.org
#
#    will be required.  Install them in advance.
#
# 2. Copy and modify example data and configuration files in the working
#    directory.
#
#    Data file: A comma- or tab-separated data sheet must have six
#       columns of 207Pb/235U, 207Pb/235U_error, 206Pb/238U,
#       206Pb/238U_error, 207Pb/206Pb, and 207Pb/206Pb_error
#
#    Configuration file: A text file describing variables used by this
#       script.  Th file name is assumed to be the same with that of
#       the data file except the extension which is cfg.
#
# [Usage]:
#
# 1. Source code:
#
#    After installation of required libraries, type like this in the
#    terminal window.
#
#    $ python UPbplot.py
#    $ python UPbplot.py -n -i data.csv
#    $ python UPbplot.py -n -i data.csv -d pdf -f
#
#    The script assumes the configuration file names is "data.cfg" as
#    defaults, if the name of the data file is "data.csv".
#
#    Command-line options:
#
#      Options:
#        -h, --help            show this help message and exit
#        -i FILE, --in=FILE    Name of input data file
#        -c FILE, --cfg=FILE   Name of configuration file
#        -o FILE, --out=FILE   Name of output file (when pdf driver is used)
#        -d DRIVER, --driver=DRIVER
#                              Choose from [pdf (default), qt5agg]
#        -f, --force-overwrite
#                              Force overwrite the pre-existing pdf
#

import os
import sys
import math
import numpy as np
import matplotlib as mpl
import pandas as pd
from json import loads
from optparse import OptionParser

# from ConfigParser import SafeConfigParser  # Python2
from configparser import ConfigParser  # Python3
from scipy import stats
from scipy import optimize
from matplotlib.patches import Ellipse

# Intersections between concordia line and error ellipses
from shapely.geometry import Polygon, LineString, LinearRing

debug = 1

# ################################################
# Initial coefficient

l238U = 1.55125 * 10 ** (-10)  # lambda_238U
l235U = 9.8485 * 10 ** (-10)  # lambda_235U
l232Th = 4.9475 * 10 ** (-11)  # lambda_232Th
l231Pa = 2.13 * 10 ** (-5)  # Rempfer2017epsl
l230Th = 9.1705 * 10 ** (-6)  # Cheng2013epsl
U85r = 137.818  # 238U/235U

# common-Pb correction
# Two-stage model for 207Pb-correction
# Stacey and Kramers (1975) EPSL
t2nd = 3.7 * 10 ** 9  # 2nd-stage startting  (3.7 Ga)
Pb64i = 11.152  # (206Pb/204Pb) at 3.7 Ga
Pb74i = 12.998  # (207Pb/204Pb) at 3.7 Ga
U8Pb4 = 9.74  # 238U/204Pb

# Time (year)
time_ka = np.array(list(range(1000, 5 * 10 ** 6, 1 * 10 ** 3)))  # 1-5000 ka
# 0.1-4601 Ma
time_ma = np.array(list(range(10 ** 5, 4600 * 10 ** 6, 10 ** 5)))

# ################################################
# Setting of file names


# ------------------------------------------------
# input filename by command line
def set_filename_input(*inputfile):
    if inputfile:
        infile = str(inputfile[0])
        if not os.path.exists(infile):
            sys.exit("Input data file %s is not found.") % infile
    else:
        print("# List of csv files")
        listfiles = os.listdir(os.getcwd())
        for i in range(len(listfiles)):
            if "csv" in listfiles[i]:
                print(listfiles[i])
        print("-------------------------")
        # infile = raw_input("Enter data filename: ")   # python2
        infile = input("Enter data filename: ")  # python3
        if not os.path.exists(infile):
            sys.exit("Input data file %s is not found.") % infile
    return infile


def set_filename_conf(filename):
    if ".cfg" in filename:
        if os.path.exists(filename):
            conffile = filename
        else:
            sys.exit("Config. data file %s is not found.") % conffile
    else:
        in_name, in_ext = os.path.splitext(filename)
        conffile = filename.replace(in_ext, ".cfg")
        if not os.path.exists(conffile):
            listfiles = os.listdir(os.getcwd())
            for i in range(len(listfiles)):
                if "cfg" in listfiles[i]:
                    print(listfiles[i])
                    print("Please enter configuration file name (*.cfg):")
                    # conffile = raw_input()  # python2
                    conffile = input()  # python3
                    if not os.path.exists(conffile):
                        sys.exit("Configuration file %s is not found.") % conffile
    return conffile


def set_filename_output(filename, driver, opt_force_overwrite):
    in_name, in_ext = os.path.splitext(filename)
    outfile = filename.replace(in_ext, ".pdf")
    if os.path.exists(outfile):
        if "pdf" in driver:
            # print(('Output file %s exists.') % outfile)
            if not (opt_force_overwrite):
                # answer = raw_input('Do you set a new file name?: [y/N] ')  #
                answer = input("Do you set a new file name?: [y/N] ")  # python3
                if (len(answer) != 0) and (answer or answer[0].lower()) == "y":
                    print("Please enter output file name (*.pdf): ")
                    # outfile = raw_input()  # python2
                    outfile = input()  # python3

    return outfile


# ################################################
# Functions

# ------------------------------------------------
# Conventional concordia diagram
# 207Pb*/235U--206Pb*/238U
def ConcLineConv(t):
    if opt_correct_disequilibrium:
        # Sakata et al., 2017, Quaternary Geochronology, eqs. (8) and (9)
        # http://dx.doi.org/10.1016/j.quageo.2016.11.002
        X = (
            np.exp(l235U * t)
            - 1
            + l235U
            / l231Pa
            * (f_Pa_U - 1)
            * (1 - np.exp(-l231Pa * t))
            * np.exp(l235U * t)
        )
        Y = (
            np.exp(l238U * t)
            - 1
            + l238U
            / l230Th
            * (f_Th_U - 1)
            * (1 - np.exp(-l230Th * t))
            * np.exp(l238U * t)
        )
    else:
        X = np.exp(l235U * t) - 1
        Y = (X + 1) ** (l238U / l235U) - 1
    return (X, Y)


# Plot a conventional concordia curve
def PlotConcConv(axcc, Xconv, Yconv, time, age_unit, L, legend_font_size):
    li = graph_label_interval * age_unit
    axcc.plot(Xconv, Yconv, color="k", linewidth=1)
    t0 = math.floor(time[0] / age_unit) * age_unit
    tann = np.arange(t0, time[-1], li)
    xa, ya = ConcLineConv(tann)

    for i, x in enumerate(xa):
        axcc.plot(xa[i], ya[i], "o", mfc="white", mec="grey", mew=1, ms=2)
        if li < age_unit:
            tlabel = "{:.1f}".format(tann[i] / age_unit)
        else:
            tlabel = "{:}".format(tann[i] / age_unit)
        axcc.annotate(
            tlabel,
            xy=(xa[i], ya[i]),
            xytext=(3, -6),
            fontsize=legend_font_size - 2,
            textcoords="offset points",
        )


# ------------------------------------------------
# Terra-Wasserburg concordia diagram
# 238U/206Pb*--207Pb*/206*Pb


def ConcLineTW(t):
    if opt_correct_disequilibrium:
        Xr = (
            np.exp(l238U * t)
            - 1
            + l238U
            / l230Th
            * (f_Th_U - 1)
            * (1 - np.exp(-l230Th * t))
            * np.exp(l238U * t)
        )
        X = 1.0 / Xr
        Y = ((1.0 / X + 1.0) ** (l235U / l238U) - 1.0) * X / U85r
    else:
        X = 1.0 / (np.exp(l238U * t) - 1)
        Y = ((1.0 / X + 1.0) ** (l235U / l238U) - 1.0) * X / U85r
    return (X, Y)


# Plot a Terra-Wasserburg concordia curve
def PlotConcTW(axtw, Xtw, Ytw, time, age_unit, L, legend_font_size):
    li = graph_label_interval * age_unit
    axtw.plot(Xtw, Ytw, color="k", linewidth=1)
    t0 = math.floor(time[0] / age_unit) * age_unit
    tann = np.arange(t0, time[-1], li)
    xa, ya = ConcLineTW(tann)
    for i, x in enumerate(xa):
        axtw.plot(xa[i], ya[i], "o", mfc="white", mec="grey", mew=1, ms=2)
        if li < age_unit:
            tlabel = "{:.1f}".format(tann[i] / age_unit)
        else:
            tlabel = "{:}".format(tann[i] / age_unit)
        axtw.annotate(
            tlabel,
            xy=(xa[i], ya[i]),
            xytext=(3, -6),
            fontsize=legend_font_size - 2,
            textcoords="offset points",
        )


def TimeRangeConv(rX, rY):
    tXmin = np.log(rX[0] + 1) / l235U
    tXmax = np.log(rX[1] + 1) / l235U
    xYmin = np.power(rY[0] + 1, 1 / (l238U / l235U)) - 1
    xYmax = np.power(rY[1] + 1, 1 / (l238U / l235U)) - 1
    tYmin = np.log(xYmin + 1) / l235U
    tYmax = np.log(xYmax + 1) / l235U
    return (tXmin, tXmax, tYmin, tYmax)


def TimeRangeTW(rx):
    txmax = np.log(1.0 / rx[0] + 1) / l238U
    txmin = np.log(1.0 / rx[1] + 1) / l238U
    return (txmin, txmax)


# ------------------------------------------------
# Age calculation of 207Pb/206Pb


def func_t76(t, j):
    res = abs(U85r * j - (np.exp(l235U * t) - 1) / (np.exp(l238U * t) - 1))
    return res


def calc_t76(age_unit, j, je, conf):
    t76 = np.empty(len(j))  # 7Pb/6Pb age
    t76upper = np.empty(len(j))  # upper 7Pb/6Pb age with error
    t76lower = np.empty(len(j))  # lower 7Pb/6Pb age
    t76_se_plus = np.empty(len(j))
    t76_se_minus = np.empty(len(j))
    cr = stats.norm.ppf(conf + (1 - conf) / 2.0)
    for i in range(len(j)):
        t76[i] = optimize.leastsq(func_t76, age_unit, args=(j[i]))[0][0]
        t76upper[i] = optimize.leastsq(func_t76, age_unit, args=(j[i] + je[i] * cr))[0][0]
        t76lower[i] = optimize.leastsq(func_t76, age_unit, args=(j[i] - je[i] * cr))[0][0]
        t76_se_plus[i] = t76upper[i] - t76[i]
        t76_se_minus[i] = t76[i] - t76lower[i]
        if t76[i] < 0.0:
            t76[i] = 0.0

    return (t76, t76_se_plus, t76_se_minus)


# ------------------------------------------------
# Calculate eigen values and vectors
def eigsorted(cov):
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    return vals[order], vecs[:, order]


def myEllipse(i, x, y, sigma_x, sigma_y, cov_xy, conf="none"):
    cov = ([sigma_x ** 2, cov_xy], [cov_xy, sigma_y ** 2])
    vals, vecs = eigsorted(cov)
    theta = np.degrees(np.arctan2(*vecs[:, 0][::-1]))
    if (vals[0] < 0.0) | (vals[1] < 0.0):
        print("!!! Unable to draw an error ellipse [Data %s] !!!" % str(i))
        ell = 0
    else:
        width, height = 2 * np.sqrt(stats.chi2.ppf(conf, 2)) * np.sqrt(vals)
        n = 100
        tpi = np.linspace(0, 2 * np.pi, n, endpoint=True)
        st = np.sin(tpi)
        ct = np.cos(tpi)
        ell = []
        angle = np.deg2rad(theta)
        sa = np.sin(angle)
        ca = np.cos(angle)
        # p = np.empty((n, 2))
        xelip = x + width / 2 * ca * ct - height / 2 * sa * st
        yelip = y + width / 2 * sa * ct + height / 2 * ca * st
    return (xelip, yelip)


# ------------------------------------------------
# One-dimentional weighted mean algorithm
# McLean et al., 2011, G-cubed, doi: 10.1029/2010GC003478


def oneWM(X, s1, conf):
    w = s1 ** (-2) / np.sum(s1 ** (-2))  # weight
    Twm = np.sum(w * X)  # weight mean of age
    S = np.sum((X - Twm) ** 2 / s1 ** 2)  # S
    # Mean Square of the Weighted Deviation
    MSWD = S / (len(X) - 1)
    # standard deviation of the weighted mean [eq. 66]
    sm = stats.norm.ppf(conf + (1 - conf) / 2.0) * np.sqrt(1.0 / np.sum(s1 ** (-2)))

    # # p.177 in Talyer1997book
    # w = s1**(-2)                # weight
    # Twm = np.sum(w*X)/np.sum(w)  # weight mean of age
    # sm = stats.norm.ppf(conf + (1 - conf)/2.)/np.sqrt(np.sum(w))
    # S = np.sum((X - Twm)**2/s1**2)  # S
    # MSWD = S/(len(X) - 1)

    return (Twm, sm, MSWD)


# ------------------------------------------------
# Two-dimentional weighted mean algorithm
# Equation numbers corresponds with those in Ludwig (1998).
def twoWM(Xi, Yi, sXi, sYi, rhoXYi, conf):
    # arguments:
    #     X, Y: measured data pont
    #    sX, sY: SD of X or Y
    N = len(Xi)
    covXYi = rhoXYi * sXi * sYi
    o11 = sYi ** 2 / ((sXi ** 2) * (sYi ** 2) - covXYi ** 2)  # eq(3)
    o22 = sXi ** 2 / ((sXi ** 2) * (sYi ** 2) - covXYi ** 2)  # eq(3)
    o12 = -covXYi / ((sXi ** 2) * (sYi ** 2) - covXYi ** 2)  # eq(3)
    # eq(6)
    x_bar = (
        np.sum(o22) * np.sum(Xi * o11 + Yi * o12)
        - np.sum(o12) * np.sum(Yi * o22 + Xi * o12)
    ) / (np.sum(o11) * np.sum(o22) - np.sum(o12) ** 2)
    # eq(7)
    y_bar = (
        np.sum(o11) * np.sum(Yi * o22 + Xi * o12)
        - np.sum(o12) * np.sum(Xi * o11 + Yi * o12)
    ) / (np.sum(o11) * np.sum(o22) - np.sum(o12) ** 2)
    # eq(2)
    (Ri, ri) = (Xi - x_bar, Yi - y_bar)
    # eq(4)
    S = np.sum(((Ri ** 2.0) * o11) + ((ri ** 2.0) * o22) + (2.0 * Ri * ri * o12))
    # eq(8)
    MSWD = S / (2 * N - 2)
    # eq(9)
    sigma_x_bar = np.sqrt(np.sum(o22) / (np.sum(o11) * np.sum(o22) - np.sum(o12) ** 2))
    # eq(9)
    sigma_y_bar = np.sqrt(np.sum(o11) / (np.sum(o11) * np.sum(o22) - np.sum(o12) ** 2))
    # eq(9)
    cov_xy_bar = -np.sum(o12) / (np.sum(o11) * np.sum(o22) - np.sum(o12) ** 2)
    return (x_bar, y_bar, MSWD, sigma_x_bar, sigma_y_bar, cov_xy_bar)


# ------------------------------------------------
# Concordia ages by using least square method
# Equation numbers corresponds with those in Ludwig (1998).


# Conventional concordia curve by Ludwig (1998)
def FitFuncConv(t, x, y, sigma_x, sigma_y, rho_xy):
    A = (x - (np.exp(l235U * t) - 1)) / sigma_x
    B = (y - (np.exp(l238U * t) - 1)) / sigma_y
    # eq(5)
    S = np.sum((A ** 2 + B ** 2 - 2 * A * B * rho_xy) / (1 - rho_xy ** 2))
    return S


# Terra-Wasserburg concordia curve by Ludwig (1998)
def FitFuncTW(t, x, y, sigma_x, sigma_y, rho_xy):
    A = (x - 1 / (np.exp(l238U * t) - 1)) / sigma_x
    B = (y - (1 / U85r) * (np.exp(l235U * t) - 1) / (np.exp(l238U * t) - 1)) / sigma_y
    # eq(5)
    S = np.sum((A ** 2 + B ** 2 - 2 * A * B * rho_xy) / (1 - rho_xy ** 2))
    return S


def ConcAgeConv(Xi, Yi, sigma_Xi, sigma_Yi, rhoXYi, Tinit=10.0 ** 6, conf=0.95):
    X_bar, Y_bar, MSWD_bar, sigma_X_bar, sigma_Y_bar, cov_XY_bar = twoWM(
        Xi, Yi, sigma_Xi, sigma_Yi, rhoXYi, conf
    )
    rho_XY_bar = cov_XY_bar / (sigma_X_bar * sigma_Y_bar)
    T_leastsq = optimize.leastsq(
        FitFuncConv, Tinit, args=(X_bar, Y_bar, sigma_X_bar, sigma_Y_bar, rho_XY_bar)
    )[0][0]
    # eq(3)
    oi = np.linalg.inv([[sigma_X_bar ** 2, cov_XY_bar], [cov_XY_bar, sigma_Y_bar ** 2]])
    # eq(14)
    Q235 = l235U * np.exp(l235U * T_leastsq)
    Q238 = l238U * np.exp(l238U * T_leastsq)
    # eq(13)
    QQ = (Q235 ** 2 * oi[0][0] + Q238 ** 2 * oi[1][1] + 2 * Q235 * Q238 * oi[0][1]) ** (
        -1
    )
    T_1sigma = np.sqrt(QQ)
    T_sigma = stats.norm.ppf(conf + (1 - conf) / 2.0) * T_1sigma
    S_bar = FitFuncConv(T_leastsq, X_bar, Y_bar, sigma_X_bar, sigma_Y_bar, rho_XY_bar)
    S = FitFuncConv(T_leastsq, Xi, Yi, sigma_Xi, sigma_Yi, rhoXYi)

    df_concordance = 1
    df_equivalence = 2 * len(Xi) - 2
    df_combined = df_concordance + df_equivalence
    MSWD_concordance = S_bar / df_concordance
    MSWD_equivalence = S / df_equivalence
    MSWD_combined = (S_bar + S) / df_combined
    P_value_eq = 1 - stats.chi2.cdf(S, df_equivalence)
    P_value_comb = 1 - stats.chi2.cdf(S_bar + S, df_combined)
    P_value_conc = 1 - stats.chi2.cdf(S_bar, df_concordance)

    return (
        T_leastsq,
        T_sigma,
        MSWD_concordance,
        MSWD_equivalence,
        MSWD_combined,
        P_value_conc,
        P_value_eq,
        P_value_comb,
    )


# Tera-Wasserburg concordia age
def ConcAgeTW(Xi, Yi, sigma_Xi, sigma_Yi, rhoXYi, Tinit=10.0 ** 6, conf=0.95):

    x = Xi  # eq. (21) = 1/(np.exp(l238U * t) - 1)
    y = Yi  # eq. (22) = 1/U85r * (np.exp(l235U * t) - 1)/(np.exp(l238U * t) - 1)
    sigma_x = sigma_Xi
    sigma_y = sigma_Yi
    rhoxy = rhoXYi
    # two-dimensional weighted mean
    x_bar, y_bar, mswd_bar, sigma_x_bar, sigma_y_bar, cov_xy_bar = twoWM(
        x, y, sigma_x, sigma_y, rhoxy, conf
    )
    rho_xy_bar = cov_xy_bar / (sigma_x_bar * sigma_y_bar)
    T_leastsq = optimize.leastsq(
        FitFuncTW, Tinit, args=(x_bar, y_bar, sigma_x_bar, sigma_y_bar, rho_xy_bar)
    )[0][0]
    # eq(3)
    oi = np.linalg.inv([[sigma_x_bar ** 2, cov_xy_bar], [cov_xy_bar, sigma_y_bar ** 2]])
    # modified from eq(14) using eq(10, 11)
    # A and B are derivative of x(t) and y(t), respectively.
    A = -l238U * np.exp(l238U * T_leastsq) / (np.exp(l238U * T_leastsq) - 1) ** 2
    B = (
        1
        / U85r
        * (
            l235U * np.exp(l235U * T_leastsq) * (np.exp(l238U * T_leastsq) - 1)
            - l238U * np.exp(l238U * T_leastsq) * (np.exp(l235U * T_leastsq) - 1)
        )
        / (np.exp(l238U * T_leastsq) - 1) ** 2
    )
    # eq(13)
    QQtw = (A ** 2 * oi[0][0] + B ** 2 * oi[1][1] + 2 * A * B * oi[0][1]) ** (-1)
    T_1sigma = np.sqrt(QQtw)
    T_sigma = stats.norm.ppf(conf + (1 - conf) / 2.0) * T_1sigma
    S_bar = FitFuncTW(T_leastsq, x_bar, y_bar, sigma_x_bar, sigma_y_bar, rho_xy_bar)
    S = FitFuncTW(T_leastsq, x, y, sigma_x, sigma_y, rhoxy)

    df_concordance = 1.0
    df_equivalence = 2.0 * len(Xi) - 2
    df_combined = df_concordance + df_equivalence
    MSWD_concordance = S_bar / df_concordance
    MSWD_equivalence = S / df_equivalence
    MSWD_combined = (S_bar + S) / df_combined
    P_value_eq = 1 - stats.chi2.cdf(S, df_equivalence)
    P_value_comb = 1 - stats.chi2.cdf(S_bar + S, df_combined)
    P_value_conc = 1 - stats.chi2.cdf(S_bar, df_concordance)

    return (
        T_leastsq,
        T_sigma,
        MSWD_concordance,
        MSWD_equivalence,
        MSWD_combined,
        P_value_conc,
        P_value_eq,
        P_value_comb,
    )


# ------------------------------------------------
# Least square method
def FitFunc_LS(x, a, b):
    return a + b * x


def FitFuncSI_LS(parameter, x, y):
    residual = y - FitFunc_LS(x, parameter[0], parameter[1])
    return residual


# Calculate slope (b) and intercept (b) by leastsq
def SlopeIntercept_LS(x, y):
    init_ab = [0.0, 0.0]
    result = optimize.leastsq(FitFuncSI_LS, init_ab, args=(x, y), full_output=1)
    parameter_optimal = result[0]
    return (parameter_optimal[1], parameter_optimal[0], result)


# ------------------------------------------------
# Maximum Likelihood Estimate (MLE)
# by York (1996) EPSL and Titterington (1979) Chem Geol


# calculate X_bar, Y_bar, and Z
# Number of equations corresponds with those in York (1969)
def Fit_XYZ(b, x, y, wx, wy, r):
    # eq(2)
    z = wx * wy / (b ** 2 * wy + wx - 2.0 * b * r * np.sqrt(wx * wy))
    # eq(4)
    x_bar = np.sum(z * x) / np.sum(z)
    y_bar = np.sum(z * y) / np.sum(z)
    a = y_bar - b * x_bar
    S = np.sum(z * (y - b * x - a) ** 2)
    return (x_bar, y_bar, z, S)


# Residual of slope (b)
def FitFuncSI(b, X, Y, wx, wy, r, case):
    X_bar, Y_bar, Z, S = Fit_XYZ(b, X, Y, wx, wy, r)
    (U, V) = (X - X_bar, Y - Y_bar)  # eq(4) in York1969epsl
    A = np.sum(Z ** 2 * ((U * V / wx) - (r * U ** 2) / np.sqrt(wx * wy)))
    B = np.sum(Z ** 2 * ((U ** 2 / wy) - (V ** 2) / wx))
    C = np.sum(Z ** 2 * ((U * V / wy) - (r * V ** 2) / np.sqrt(wx * wy)))
    if case == 1:
        # eq(2) in Titterington1979chemg
        S = (-B + np.sqrt(B ** 2 + 4 * A * C)) / (2 * A) - b
    else:
        S = (-B - np.sqrt(B ** 2 + 4 * A * C)) / (2 * A) - b
    return S


# Calculate slope (b) and intercept (b) by MLE
def SlopeIntercept(x, y, sigma_x, sigma_y, rho_xy, case):
    wx = sigma_x ** (-2)
    wy = sigma_y ** (-2)
    init_b = 1.0

    # compare two solutions of b
    if case == 1:
        bi1 = optimize.leastsq(FitFuncSI, init_b, args=(x, y, wx, wy, rho_xy, case))[0][0]
        X_bar1, Y_bar1, Z1, S1 = Fit_XYZ(bi1, x, y, wx, wy, rho_xy)
    else:
        bi2 = optimize.leastsq(FitFuncSI, init_b, args=(x, y, wx, wy, rho_xy, case))[0][0]
        X_bar2, Y_bar2, Z2, S2 = Fit_XYZ(bi2, x, y, wx, wy, rho_xy)

    # Compare the residuals between two solituion
    if case == 1:
        X_bar, Y_bar, Z, bi = X_bar1, Y_bar1, Z1, bi1
    else:
        X_bar, Y_bar, Z, bi = X_bar2, Y_bar2, Z2, bi2

    # In practice the positive square root will be appropriate.
    # (p. 185 in Titterngton1979chemg)
    # X_bar, Y_bar, Z, bi = X_bar1, Y_bar1, Z1, bi1

    ai = Y_bar - bi * X_bar

    # Appendix in Titterngton1979chemg
    sigma_bi = np.sqrt(1.0 / np.sum(Z * (x - X_bar) ** 2))
    # sigma_ai = np.sqrt(np.sum(x**2*Z)/(np.sum(Z)*np.sum(Z*(x-X_bar)**2)))
    sigma_ai = np.sqrt(1.0 / np.sum(Z) + X_bar ** 2 * sigma_bi ** 2)  # York2004ajp
    return (X_bar, Y_bar, ai, bi, sigma_ai, sigma_bi)


# ------------------------------------------------
# Intercept age on conventional concordia diagram
# Ludwig (1980) EPSL


# 1 sigma of straight line
def SIsigma(x, x_bar, y_bar, b, sigma_a, sigma_b, conf=0.95):
    # eq(29)
    sigma_a2 = stats.norm.ppf(conf + (1 - conf) / 2.0) * sigma_a
    sigma_b2 = stats.norm.ppf(conf + (1 - conf) / 2.0) * sigma_b
    sigma2 = np.sqrt(sigma_a2 ** 2.0 + sigma_b2 ** 2.0 * x * (x - 2.0 * x_bar))
    return sigma2


def SIsigma2(x, x_bar, y_bar, b, sigma_a, sigma_b, conf=0.95):
    # eq(30)
    dtheta = sigma_b * np.cos(np.arctan(b)) ** 2.0
    b2 = (np.tan(np.arctan(b) + dtheta) + np.tan(np.arctan(b) - dtheta)) / 2.0
    a2 = y_bar - b2 * x_bar
    sigma_b2 = (np.tan(np.arctan(b) + dtheta) - np.tan(np.arctan(b) - dtheta)) / 2.0
    sigma_a2 = sigma_a + (sigma_b2 - sigma_b) * x_bar
    sigma_b2 = stats.norm.ppf(conf + (1 - conf) / 2.0) * sigma_b2
    sigma_a2 = stats.norm.ppf(conf + (1 - conf) / 2.0) * sigma_a2
    # eq(31)
    sigma2 = np.sqrt(sigma_a2 ** 2 + sigma_b2 ** 2.0 * x * (x - 2.0 * x_bar))
    return (a2, sigma_a2, b2, sigma_b2, sigma2)


# Intercept age on conventional diagram
def FitFuncSIageConv(t, a, b, sigma_a, sigma_b, x_bar, y_bar, conf, case):
    x1 = np.exp(l235U * t) - 1.0
    y1 = np.exp(l238U * t) - 1.0
    sigma = SIsigma(x1, x_bar, y_bar, b, sigma_a, sigma_b, conf=conf)
    if case == 1:
        y2 = b * x1 + a - sigma
    elif case == 2:
        y2 = b * x1 + a + sigma
    else:
        y2 = b * x1 + a
    S = (y2 - y1) ** 2
    return S


def SIageConv(a, b, sigma_a, sigma_b, x_bar, y_bar, init_t=10 ** 6, conf=0.95):
    T = optimize.leastsq(
        FitFuncSIageConv, init_t, args=(a, b, sigma_a, sigma_b, x_bar, y_bar, conf, 0)
    )[0][0]
    Tmin = optimize.leastsq(
        FitFuncSIageConv, init_t, args=(a, b, sigma_a, sigma_b, x_bar, y_bar, conf, 1)
    )[0][0]
    Tmax = optimize.leastsq(
        FitFuncSIageConv, init_t, args=(a, b, sigma_a, sigma_b, x_bar, y_bar, conf, 2)
    )[0][0]
    if Tmax < Tmin:
        Tmax, Tmin = Tmin, Tmax
    return (T, Tmin, Tmax)


# ------------------------------------------------
# Intercept age on Tera-Wasserburg diagram
def FitFuncSIageTW(t, a, b, sigma_a, sigma_b, x_bar, y_bar, conf, cs):
    x1 = 1.0 / (np.exp(l238U * t) - 1)
    y1 = 1.0 / U85r * (np.exp(l235U * t) - 1) / (np.exp(l238U * t) - 1)
    sigma = SIsigma(x1, x_bar, y_bar, b, sigma_a, sigma_b, conf=conf)
    if cs == 1:
        y2 = b * x1 + a - sigma
    elif cs == 2:
        y2 = b * x1 + a + sigma
    else:
        y2 = b * x1 + a
    S = (y2 - y1) ** 2
    return S


def SIageTW(a, b, sigma_a, sigma_b, x_bar, y_bar, init_t=10 ** 6, conf=0.95):
    T = optimize.leastsq(
        FitFuncSIageTW, init_t, args=(a, b, sigma_a, sigma_b, x_bar, y_bar, conf, 0)
    )[0][0]
    Tmin = optimize.leastsq(
        FitFuncSIageTW, init_t, args=(a, b, sigma_a, sigma_b, x_bar, y_bar, conf, 1)
    )[0][0]
    Tmax = optimize.leastsq(
        FitFuncSIageTW, init_t, args=(a, b, sigma_a, sigma_b, x_bar, y_bar, conf, 2)
    )[0][0]
    if Tmax < Tmin:
        Tmax, Tmin = Tmin, Tmax
    return (T, Tmin, Tmax)


# ------------------------------------------------
# Weighted mean and reduced chi-square
# Spencer2016gf
def calc_chi2_red(x, s1, wm, n, opt):
    # x: ages of each sample
    # s1: errors (1 sigma) of each sample
    # wm: weighted mean age
    # n: number of sample
    chi2_red = 1 / (n - 1) * np.sum((x - wm) ** 2 / s1 ** 2)
    error_min = 1 - 2 * np.sqrt(2 / (n - 1))
    error_max = 1 + 2 * np.sqrt(2 / (n - 1))
    if (chi2_red <= error_max) & (chi2_red >= error_min):
        res = "Passed"
    else:
        if opt:
            if chi2_red < error_min:
                res = "Failed, under dispersion or uncertainties overestimated"
            else:
                res = "Failed, over dispersion or uncertainties underestimated"
        else:
            res = "Failed"
    return (chi2_red, res)


# ------------------------------------------------
# Exclude outlier: Generalized ESD test: Rosner 1983
# ESD (extreme Studentized deviate)
#
# https://www.itl.nist.gov/div898/handbook/eda/section3/eda35h3.htm
#
# Rosner, Bernard (May 1983), Percentage Points for a Generalized ESD
#            Many-Outlier Procedure,Technometrics, 25(2), pp. 165-172.
def GESDtest(Tall, s1, ind, cr):

    conf = 1 - cr  # significant level is conf (e.g. 0.05)
    ss = 2.0  # 2 sigma level
    x0 = Tall[ind]
    x = x0
    s0 = s1[ind]
    s = s0
    ii = ind
    oo = []
    while True:
        x = Tall[ii]
        s = s1[ii]
        n = len(x)
        t = stats.t.ppf(q=(1 - cr / (2 * n)), df=n - 2)
        tau = (n - 1) * t / np.sqrt(n * (n - 2) + n * t * t)
        xs_min, xs_max = np.min(x + ss * s), np.max(x - ss * s)

        i_max = [i for i in ind if x0[i] - ss * s0[i] == xs_max]
        # if data have the multiple maximum values, the datum with large error will be excluded.
        if len(i_max) > 1:
            smax = s0.min()
            for j in i_max:
                if s0[j] > smax:
                    smax = s0[j]
            ii_max = [i for i in ind if s0[i] == smax]

        i_min = [i for i in ind if x0[i] + ss * s0[i] == xs_min]
        # if data have the multiple minimum values, the datum with larger error will be excluded.
        if len(i_min) > 1:
            smax = s0.min()
            for j in i_min:
                if s0[j] > smax:
                    smax = s0[j]
            i_min = [i for i in ind if s0[i] == smax]

        Twm, sm, MSWD = oneWM(x[ii], s[ii], conf)
        if np.abs(xs_max - Twm) > np.abs(xs_min - Twm) and xs_max >= Twm:
            i_far = i_max[-1]
            tau_far = np.abs(((x[i_far] - ss * s[i_far]) - Twm) / sm)
        elif np.abs(xs_max - Twm) <= np.abs(xs_min - Twm) and xs_min <= Twm:
            i_far = i_min[-1]
            tau_far = np.abs(((x[i_far] + ss * s[i_far]) - Twm) / sm)
        else:
            break

        if float(tau_far) < tau:
            break

        ii = np.delete(ii, np.where(ii == i_far))

    oo = np.setdiff1d(ind, ii)
    return (ii, oo)


# ------------------------------------------------
# range calculation
def calc_legend_pos(range_XY):
    legend_pos_x = [0.02] * 10
    legend_pos_y = [0] * 10
    for i in range(0, 10):
        legend_pos_y[i] = 0.98 - 0.05 * i

    return (legend_pos_x, legend_pos_y)


# ------------------------------------------------
# Legend
def legend_data_number(ax, axn, x, y, ind):
    # ax[axn].text(x, y, '$N$ = %d, $n$ = %d' % (N, n_in),
    #             fontsize=legend_font_size)
    n_in = len(ind)
    ax[axn].text(
        0.02,
        0.98,
        "$N$ = %d, $n$ = %d" % (N, n_in),
        transform=ax[axn].transAxes,
        verticalalignment="top",
        fontsize=legend_font_size,
    )


# ------------------------------------------------
# plot data points
def plot_data_point(ax, X, Y, ind, outd, outd_disc):
    for i in range(len(X)):
        if i in ind:
            ax.plot(
                X[i],
                Y[i],
                dp1_marker_type,
                alpha=dp1_marker_alpha,
                markersize=dp1_marker_size,
                markerfacecolor=dp1_marker_fc,
                markeredgecolor=dp1_marker_ec,
                markeredgewidth=dp1_marker_ew,
            )
        elif i in outd_disc:
            ax.plot(
                X[i],
                Y[i],
                dp1_marker_type,
                alpha=dp2_marker_alpha,
                markersize=dp2_marker_size,
                markerfacecolor=dp2_marker_fc,
                markeredgecolor=dp2_marker_ec,
                markeredgewidth=dp2_marker_ew,
            )
        else:
            ax.plot(
                X[i],
                Y[i],
                dp0_marker_type,
                alpha=dp0_marker_alpha,
                markersize=dp0_marker_size,
                markerfacecolor=dp0_marker_fc,
                markeredgecolor=dp0_marker_ec,
                markeredgewidth=dp0_marker_ew,
            )


# ------------------------------------------------
# Intersections between concordia line and error ellipses
def intersections_ellipse(x, y, line):
    # Intersections between concordia line and error ellipses
    ellipse = []
    yl = []
    for i, j in enumerate(x):
        ellipse.append([x[i], y[i]])

    shapely_poly = LinearRing(ellipse)
    shapely_line = LineString(line)
    mp = shapely_poly.intersection(shapely_line)

    if mp.is_empty:
        return 0
    elif mp.geom_type == "Point":
        return 1
    elif mp.geom_type == "MultiPoint":
        return 1
    else:
        return 0


# ------------------------------------------------
# Judgement of discordance (disc_type == 5)
def discordant_judge(xd, yd, sigma_xd, sigma_yd, cov_xdyd, conf, i_in, i_disc, line):
    i_ind = i_in
    for j, i in enumerate(i_in):
        dp_ell_x, dp_ell_y = myEllipse(
            i, xd[i], yd[i], sigma_xd[i], sigma_yd[i], cov_xdyd[i], conf,
        )
        if not intersections_ellipse(dp_ell_x, dp_ell_y, line):
            i_disc = np.append(i_disc, i)

    i_disc = np.unique(i_disc)
    i_ind = np.setdiff1d(i_in, i_disc)
    return (i_ind, i_disc)


# ------------------------------------------------
# Discordance (%) = [1 - A/B] * 100
def discordance(
    t75, t75e, t68, t68e, t76, t76_min, t76_max, sd, method,
):

    if method == 0:
        disc_percent = (1 - (t68 / t76)) * 100.0
    elif method == 1:
        disc_percent = (1 - (t75 / t76)) * 100.0
    elif method == 2:
        disc_percent = (1 - (t68 / t75)) * 100.0
    elif method == 3:
        disc_percent = (1 - (t75 / t68)) * 100.0
    elif method == 4:
        disc_percent = (1 - ((t75 - t75e * sd) / (t68 + t68e * sd))) * 100.0
    else:
        sys.exit("Exit at function of discordance")

    return disc_percent


# ------------------------------------------------
# print discordant data
def print_discordant_data(tdisc, odisc, pdisc, sigma):
    # tdisc: disc_type
    # odisc: outd_disc
    # pdisc: disc_percent
    if tdisc == 0:
        print("Discordance is calculated by", end=" ")
        print("100*(1-([206Pb/238U age]/[207Pb/206Pb age]))")
        # Discordant data points
        print("Discordant data points [n = %d] are" % len(odisc))
        for i in odisc:
            print(
                "%d: %s%% = (1-%.1f/%.1f) x 100"
                % (
                    i + 1,
                    format(pdisc[i], dignum),
                    d["t68"].iloc[i] / age_unit,
                    d["t75"].iloc[i] / age_unit,
                )
            )
    elif tdisc == 1:
        print("Discordance is calculated by", end=" ")
        print("100*(1-([207Pb/235U age]/[207Pb/206Pb age]))")
        # Discordant data points
        print("Discordant data points [n = %d] are" % len(odisc))
        for i in odisc:
            print(
                "%d: %s%% = (1-%.1f/%.1f) x 100"
                % (
                    i + 1,
                    format(pdisc[i], dignum),
                    d["t75"].iloc[i] / age_unit,
                    d["t76"].iloc[i] / age_unit,
                )
            )
    elif tdisc == 2:
        print("Discordance is calculated by", end=" ")
        print("100*(1-([206Pb/238U age]/[207Pb/235U age])")
        # Discordant data points
        print("Discordant data points [n = %d] are" % len(odisc))
        for i in odisc:
            print(
                "%d: %s%% = (1-%.1f/%.1f) x 100"
                % (
                    i + 1,
                    format(pdisc[i], dignum),
                    d["t68"].iloc[i] / age_unit,
                    d["t75"].iloc[i] / age_unit,
                )
            )
    elif tdisc == 3:
        print("Discordance is calculated by", end=" ")
        print("100*(1-([207Pb/235U age]/[206Pb/238U age])")
        # Discordant data points
        print("Discordant data points [n = %d] are" % len(odisc))
        for i in odisc:
            print(
                "%d: %s%% = (1-%.1f/%.1f) x 100"
                % (
                    i + 1,
                    format(pdisc[i], dignum),
                    d["t75"].iloc[i] / age_unit,
                    d["t68"].iloc[i] / age_unit,
                )
            )
    elif tdisc == 4:
        print("Discordance is calculated by", end=" ")
        print("100*(1-(min[207Pb/235U age] / max[206Pb/238U age])")
        # Discordant data points
        print("Discordant data points [n = %d] are" % len(odisc))
        for i in odisc:
            print(
                "%d: %s%% = (1-%.1f/%.1f) x 100"
                % (
                    i + 1,
                    format(pdisc[i], dignum),
                    (d["t75"].iloc[i] - d["t75e"].iloc[i] * sigma) / age_unit,
                    (d["t68"].iloc[i] + d["t68e"].iloc[i] * sigma) / age_unit,
                )
            )


# ------------------------------------------------
# ages for correction for initial disequilibria
# eqs. (1) and (2) in Sakata (2018) Geochemical Journal
def func_Tdiseq(t, r, rtype):
    if rtype == 68:
        F = abs(
            (np.exp(l238U * t) - 1)
            + l238U
            / l230Th
            * (f_Th_U - 1)
            * (1 - np.exp(-l230Th * t))
            * np.exp(l238U * t)
            - r
        )
    elif rtype == 75:
        F = abs(
            (np.exp(l235U * t) - 1)
            + l235U
            / l231Pa
            * (f_Pa_U - 1)
            * (1 - np.exp(-l231Pa * t))
            * np.exp(l235U * t)
            - r
        )
    return F


def SI_Tdiseq(R75m, R68m, rtype):
    T = []  # corrected age for initial disequilibria
    if rtype == 68:
        for i, r in enumerate(R68m):
            t = 1 / l238U * np.log(r + 1)
            T68 = optimize.leastsq(func_Tdiseq, t, args=(r, rtype))[0][0]
            T.append(T68)
    elif rtype == 75:
        for i, r in enumerate(R75m):
            t = 1 / l235U * np.log(r + 1)
            T75 = optimize.leastsq(func_Tdiseq, t, args=(r, rtype))[0][0]
            T.append(T75)
    else:
        print("Please define rtype = 68 or 75.")
        exit
    return T


# ------------------------------------------------
# Common-Pb correction: 207Pb method

# Zricon crystalization age from eq (5) in Sakata2018GeochemJ
def func_tPb76c(t, R76c, R76m, R86m):
    F = abs(
        (
            1
            / U85r
            * (
                (np.exp(l235U * t) - 1)
                + l235U
                / l231Pa
                * (f_Pa_U - 1)
                * (1 - np.exp(-l235U * t))
                * np.exp(l235U * t)
            )
            - (R76m - R76c) / R86m
        )
        / (
            (np.exp(l238U * t) - 1)
            + l238U
            / l230Th
            * (f_Th_U - 1)
            * (1 - np.exp(-l230Th * t))
            * np.exp(-l238U * t)
        )
        - R76c
    )
    return F


def SI_Pb76c(R76c, R76m, R86m):
    t = 1 / l238U * np.log(1 / R86m + 1)
    T = optimize.leastsq(func_tPb76c, t, args=(R76c, R76m, R86m))[0][0]
    return T


def func_Pb76c(X, Y):
    "Calculate common 207Pb/206Pb at t1 age"
    # y0 = common 207Pb/206Pb at time t1
    # x1 = measured 238U/206Pb
    # y1 = measured 207Pb/206Pb
    # x2 = corrected (radiogenic) 238U/206Pb
    # y2 = corrected (radiogenic) 207Pb/206Pb
    xx2 = []  # corrected (radiogenic) 238U/206Pb
    yy2 = []  # corrected (radiogenic) 207Pb/206Pb
    f206p = []  # fraction of common 206Pb
    T = []  # zircon crystalization age from eq. (5) in Sakata2018GeochemJ
    for i, y1 in enumerate(Y):
        x1 = X[i]
        t1 = 1 / l238U * np.log(1 / x1 + 1)
        Pb64c_t1 = Pb64i + U8Pb4 * (np.exp(l238U * t2nd) - np.exp(l238U * t1))
        Pb74c_t1 = Pb74i + U8Pb4 / U85r * (np.exp(l235U * t2nd) - np.exp(l235U * t1))
        Pb76c0 = Pb74c_t1 / Pb64c_t1
        x0 = 0
        y0 = Pb76c0
        line0 = LineString([(0.0, y0), (10 ** 5, y0)])
        Xtw, Ytw = ConcLineTW(np.array(list(range(10 ** 5, 5 * 10 ** 9, 10 ** 5))))
        line_tw = LineString([(Xtw[i], Ytw[i]) for i, j in enumerate(Xtw)])
        res = line0.intersection(line_tw)
        # y = a x + b
        # a = delta(207Pb/206Pb)/delta(238U/206Pb)
        a = (y1 - y0) / (x1 - x0)
        b = y1 - a * x1
        x3 = -b / a
        line1 = LineString([[x0, y0], [x1, y1], [x3, 0]])
        res1 = line1.intersection(line_tw)
        x2 = res1[1].x  # np.max(res1.bounds)
        y2 = res1[1].y  # np.min(res1.bounds)
        f = (y1 - y2) / (y0 - y2) * 100
        t68c = SI_Pb76c(y0, y1, x1)  # y0 = R76c, y1 = R76m, x1 = R86m
        t68 = 1 / l238U * np.log(1 / x1 + 1)

        xx2.append(x2)
        yy2.append(y2)
        f206p.append(f)
        if opt_correct_disequilibrium:
            T.append(t68c)
        else:
            T.append(t68)

    return (xx2, yy2, f206p, T)


# ------------------------------------------------
# draw error ellipses of data points
def plot_data_point_error_ellipse_disc(
    ax, axn, X, Y, sigma_X, sigma_Y, cov_XY, cr, ind, outd, outd_disc, line
):
    for i in range(len(X)):
        dp_ell_x, dp_ell_y = myEllipse(
            i, X[i], Y[i], sigma_X[i], sigma_Y[i], cov_XY[i], conf=cr,
        )

        if i in outd:
            ax[axn].plot(
                dp_ell_x,
                dp_ell_y,
                alpha=dp0_ee_alpha,
                linestyle=dp0_ee_ls,
                color=dp0_ee_ec,
                linewidth=dp0_ee_ew,
            )
        elif i in ind:
            ax[axn].plot(
                dp_ell_x,
                dp_ell_y,
                alpha=dp1_ee_alpha,
                linestyle=dp1_ee_ls,
                color=dp1_ee_ec,
                linewidth=dp1_ee_ew,
            )
        elif i in outd_disc:
            ax[axn].plot(
                dp_ell_x,
                dp_ell_y,
                alpha=dp2_ee_alpha,
                linestyle=dp2_ee_ls,
                color=dp2_ee_ec,
                linewidth=dp2_ee_ew,
            )

    return (ind, outd_disc)


# ------------------------------------------------
# draw error ellipses of data points
def plot_data_point_error_ellipse(
    ax, axn, X, Y, sigma_X, sigma_Y, cov_XY, cr, ind, outd, outd_disc
):
    for i in range(len(X)):
        dp_ell_x, dp_ell_y = myEllipse(
            i, X[i], Y[i], sigma_X[i], sigma_Y[i], cov_XY[i], conf=cr,
        )
        if i in ind:
            ax[axn].plot(
                dp_ell_x,
                dp_ell_y,
                alpha=dp1_ee_alpha,
                linestyle=dp1_ee_ls,
                color=dp1_ee_ec,
                linewidth=dp1_ee_ew,
            )
        elif i in outd_disc:
            ax[axn].plot(
                dp_ell_x,
                dp_ell_y,
                alpha=dp2_ee_alpha,
                linestyle=dp2_ee_ls,
                color=dp2_ee_ec,
                linewidth=dp2_ee_ew,
            )
        else:
            ax[axn].plot(
                dp_ell_x,
                dp_ell_y,
                alpha=dp0_ee_alpha,
                linestyle=dp0_ee_ls,
                color=dp0_ee_ec,
                linewidth=dp0_ee_ew,
            )


# ------------------------------------------------
# Plot two-dimensional weighted mean ages and confidence ellipses
def plot_2D_wm(ax, axn, X, Y, sigma_X, sigma_Y, rho_XY, cr, legend_pos_x, legend_pos_y):
    Xwm_bar, Ywm_bar, MSWDwm, sigma_Xwm_bar, sigma_Ywm_bar, cov_XYwm_bar = twoWM(
        X, Y, sigma_X, sigma_Y, rho_XY, conf=cr
    )
    twm_ell_x, twm_ell_y = myEllipse(
        0, Xwm_bar, Ywm_bar, sigma_Xwm_bar, sigma_Ywm_bar, cov_XYwm_bar, conf=cr,
    )
    ax[axn].fill(
        twm_ell_x, twm_ell_y, fc=twm_ee_fc, alpha=twm_ee_alpha,
    )
    ax[axn].plot(
        twm_ell_x, twm_ell_y, color=twm_ee_ec, linewidth=twm_ee_ew,
    )

    # legend
    ax[axn].text(
        legend_pos_x,
        legend_pos_y,
        "2D weighted mean [%d%% conf.] (MSWD=%s)" % (cr * 100, format(MSWDwm, dignum)),
        transform=ax[axn].transAxes,
        verticalalignment="top",
        fontsize=legend_font_size,
    )


# ------------------------------------------------
# concordia ages on concordia lines
def concordia_age(ctype, X, Y, sigma_X, sigma_Y, rho_XY, cr):
    if ctype == "conv":
        T_lsq, S_lsq, MSWDconc, MSWDeq, MSWDcomb, Pconc, Peq, Pcomb = ConcAgeConv(
            X, Y, sigma_X, sigma_Y, rho_XY, Tinit=age_unit, conf=cr
        )
        X_lsq, Y_lsq = ConcLineConv(T_lsq)
    elif ctype == "tw":
        T_lsq, S_lsq, MSWDconc, MSWDeq, MSWDcomb, Pconc, Peq, Pcomb = ConcAgeTW(
            X, Y, sigma_X, sigma_Y, rho_XY, Tinit=age_unit, conf=cr
        )
        X_lsq, Y_lsq = ConcLineTW(T_lsq)

    else:
        sys.exit('Please choose type of concordia, "conv" or "tw"')
    return (T_lsq, S_lsq, X_lsq, Y_lsq, MSWDconc, MSWDeq, MSWDcomb, Pconc, Peq, Pcomb)


# ------------------------------------------------
# Plot concordia ages on concordia lines
def plot_concordia_age(
    ax, axn, T_lsq, S_lsq, X_lsq, Y_lsq, MSWD, cr, legend_pos_x, legend_pos_y
):

    ax[axn].plot(
        X_lsq,
        Y_lsq,
        ca_marker_type,
        color=ca_marker_fc,
        markersize=ca_marker_size,
        markeredgecolor=ca_marker_ec,
        markeredgewidth=1.0,
    )

    # legend
    ax[axn].text(
        legend_pos_x,
        legend_pos_y,
        u"Concordia age = %s ± %s %s [%d%% conf.]"
        % (
            format(T_lsq / age_unit, dignum),
            format(S_lsq / age_unit, dignum),
            age_unit_name,
            cr * 100,
        ),
        transform=ax[axn].transAxes,
        verticalalignment="top",
        fontsize=legend_font_size,
    )


# ------------------------------------------------
# Text for MSWD of concordia ages
def plot_concordia_age_MSWD(ax, axn, MSWD, ca_mswd, p_value, legend_pos_x, legend_pos_y):

    if ca_mswd == 1:
        ax[axn].text(
            legend_pos_x,
            legend_pos_y,
            r"(MSWD of equivalence = %s, p($\chi^2$) = %s)"
            % (format(MSWD, dignum), format(p_value, dignum2)),
            transform=ax[axn].transAxes,
            verticalalignment="top",
            fontsize=legend_font_size,
        )
    elif ca_mswd == 2:
        ax[axn].text(
            legend_pos_x,
            legend_pos_y,
            r"(MSWD of combined = %s, p($\chi^2$) = %s)"
            % (format(MSWD, dignum), format(p_value, dignum2)),
            transform=ax[axn].transAxes,
            verticalalignment="top",
            fontsize=legend_font_size,
        )
    else:
        ax[axn].text(
            legend_pos_x,
            legend_pos_y,
            r"(MSWD of concordance = %s, p($\chi^2$) = %s)"
            % (format(MSWD, dignum), format(p_value, dignum2)),
            transform=ax[axn].transAxes,
            verticalalignment="top",
            fontsize=legend_font_size,
        )


# ------------------------------------------------
# Plot concordia-intercept lines and confidence band
def plot_concordia_intercept_age(
    ax,
    axn,
    ctype,
    X,
    Y,
    sigma_X,
    sigma_Y,
    cr,
    rho_XY,
    range_XY,
    T_lsq,
    case,
    legend_pos_x,
    legend_pos_y,
):
    xx = np.linspace(range_XY[0][0], range_XY[0][1], 5000)
    Xsi_bar, Ysi_bar, ai, bi, sigma_a, sigma_b = SlopeIntercept(
        X, Y, sigma_X, sigma_Y, rho_XY, case
    )

    # Intercept age
    if ctype == "conv":
        Tsi, Tmin, Tmax = SIageConv(
            ai, bi, sigma_a, sigma_b, Xsi_bar, Ysi_bar, init_t=T_lsq, conf=cr
        )
    elif ctype == "tw":
        Tsi, Tmin, Tmax = SIageTW(
            ai, bi, sigma_a, sigma_b, Xsi_bar, Ysi_bar, init_t=T_lsq, conf=cr
        )

    # confidence band
    sigma = SIsigma(xx, Xsi_bar, Ysi_bar, bi, sigma_a, sigma_b, conf=cr)
    y1 = bi * xx + ai - sigma
    y2 = bi * xx + ai + sigma
    # a2, sigma_a2, b2, sigma_b2, sigma2 = SIsigma2(
    #     xx, Xsi_bar, Ysi_bar, bi, sigma_a, sigma_b, conf=cr)
    # y1 = b2*xx+a2-sigma2
    # y2 = b2*xx+a2+sigma2

    ax[axn].fill_between(
        xx,
        y1,
        y2,
        where=y2 >= y1,
        facecolor=ia_fill_color,
        edgecolor="None",
        alpha=ia_alpha,
        interpolate=True,
    )
    ax[axn].plot(xx, bi * xx + ai, linewidth=ia_line_width, color=ia_line_color)

    # legend
    ax[axn].text(
        legend_pos_x,
        legend_pos_y,
        "Intercept age = %s +%s %s %s (%d%% conf.)"
        % (
            format(Tsi / age_unit, dignum),
            format((Tmax - Tsi) / age_unit, dignum),
            format((Tmin - Tsi) / age_unit, dignum),
            age_unit_name,
            cr * 100,
        ),
        transform=ax[axn].transAxes,
        verticalalignment="top",
        fontsize=legend_font_size,
    )

    return (Tsi, Tmax, Tmin)


# ------------------------------------------------
# Choose one age type among 206Pb/238U, 207Pb/235U, or 207Pb/206Pb
def select_age_type(age_type):
    if age_type == 68:
        Tall = d["t68"] / age_unit
        s1 = Tall * SY
        label = "$^{206}$Pb* / $^{238}$U age (%s)" % age_unit_name
    elif age_type == 75:
        Tall = d["t75"] / age_unit
        s1 = Tall * SX
        label = "$^{207}$Pb* / $^{235}$U age (%s)" % age_unit_name
    elif age_type == 76:
        Tall = d["t76"] / age_unit
        s1 = Tall * Sy
        label = "$^{207}$Pb* / $^{206}$Pb* age (%s)" % age_unit_name
    else:
        sys.exit("Error at select_age_type in one-dimensional bar plot.")

    return (Tall, s1, label)


# ------------------------------------------------
# Plot one-dimensional weighted mean, SD, and MSWD
def plot_oneD_weighted_mean(
    ax_1D, oneD_age_type, Tall, s1, ind, outd, outd_disc, cr, legend_pos_x, legend_pos_y
):

    if oneD_yaxis_log == 1:
        ax_1D.set_yscale("log")

    Twm, sm, MSWD = oneWM(Tall[ind], s1[ind], conf=cr)

    # confidence band of the weighted mean
    ax_1D.axhspan(Twm - sm, Twm + sm, facecolor=oneD_band_fc, alpha=oneD_band_alpha)

    # plot 2D weighted mean
    ax_1D.plot(
        [0.0, len(Tall) + 1],
        [Twm, Twm],
        linewidth=oneD_wm_line_width,
        color=oneD_wm_line_color,
    )

    # Sorting age data
    if oneD_plot_sort == 1:
        Tplot = np.sort(Tall)
        Tnumd = np.argsort(Tall)
    else:
        Tplot = Tall
        Tnumd = np.where(Tall)[0]

    ind_sort = []
    outd_sort = []
    outd_disc_sort = []
    for n, t in enumerate(Tplot):
        i = Tnumd[n]
        if len(np.intersect1d(i, ind)) > 0:
            eb1 = ax_1D.errorbar(
                n + 1,
                t,
                yerr=stats.norm.ppf(cr + (1 - cr) / 2.0) * s1[i],
                ecolor=dp1_bar_color,
                linewidth=dp1_bar_line_width,
            )
            eb1[-1][0].set_linestyle(dp1_bar_line_style)
            ind_sort.append(n)
        elif len(np.intersect1d(i, outd_disc)) > 0:
            eb2 = ax_1D.errorbar(
                n + 1,
                t,
                yerr=stats.norm.ppf(cr + (1 - cr) / 2.0) * s1[i],
                ecolor=dp2_bar_color,
                linewidth=dp2_bar_line_width,
            )
            eb2[-1][0].set_linestyle(dp2_bar_line_style)
            outd_disc_sort.append(n)
        else:
            eb0 = ax_1D.errorbar(
                n + 1,
                t,
                yerr=stats.norm.ppf(cr + (1 - cr) / 2.0) * s1[i],
                ecolor=dp0_bar_color,
                linewidth=dp0_bar_line_width,
            )
            eb0[-1][0].set_linestyle(dp0_bar_line_style)
            outd_sort.append(n)

    # plot data point
    if opt_data_point:
        plot_data_point(
            ax_1D, np.where(Tall)[0] + 1, Tplot, ind_sort, outd_sort, outd_disc_sort
        )

    # legend
    legend_data_number(ax, axn, legend_pos_x[0], legend_pos_y[0], ind)

    oneD_xticks = np.arange(1, len(Tall) + 1, 1)
    ax_1D.set_xticks(oneD_xticks)
    ax_1D.set_xticklabels(Tnumd + 1, rotation=270, fontsize="small")

    ax_1D.text(
        legend_pos_x[0],
        legend_pos_y[1],
        u"Weighted mean = %s ± %s %s [%d%% conf.] (MSWD = %s)"
        % (
            format(Twm, dignum),
            format(sm, dignum),
            age_unit_name,
            cr * 100,
            format(MSWD, dignum),
        ),
        transform=ax_1D.transAxes,
        verticalalignment="top",
        fontsize=legend_font_size,
    )

    ax_1D.text(
        legend_pos_x[0],
        legend_pos_y[2],
        "Error bars are %d%% conf." % (cr * 100),
        transform=ax_1D.transAxes,
        verticalalignment="top",
        fontsize=legend_font_size,
    )

    chi2_red, res_chi2_red = calc_chi2_red(Tall[ind], s1[ind], Twm, len(ind), opt=0)

    ax_1D.text(
        legend_pos_x[0],
        legend_pos_y[3],
        "$\chi^2_{red}$ = %s (%s)" % (format(chi2_red, dignum), res_chi2_red),
        transform=ax[axn].transAxes,
        verticalalignment="top",
        fontsize=legend_font_size,
    )

    return (Twm, sm, MSWD)


# ------------------------------------------------
# Th/U vs age plot
# Descriminating between metamorphic and igneous origins
# Uncetainty in Th/U is assumed to be ± 10%
def plot_Th_U(
    axb, Th_U, Th_U_e, Tall, s1, ind, outd, outd_disc, cr, range_hist_x, range_hist_y2
):
    axb.set_xlim(range_hist_x[0], range_hist_x[1])
    axb.set_ylim(range_hist_y2[0], range_hist_y2[1])
    axb.set_ylabel("Th/U", fontsize=legend_font_size + 4)
    for i in range(len(Tall)):
        if i in outd:
            eb0 = axb.errorbar(
                Tall[i],
                Th_U[i],
                xerr=stats.norm.ppf(cr + (1 - cr) / 2.0) * s1[i],
                yerr=Th_U_e[i],
                ecolor=dp0_bar_color,
                linewidth=dp0_bar_line_width,
            )
            eb0[-1][0].set_linestyle(dp0_bar_line_style)
        elif i in outd_disc:
            eb2 = axb.errorbar(
                Tall[i],
                Th_U[i],
                xerr=stats.norm.ppf(cr + (1 - cr) / 2.0) * s1[i],
                yerr=Th_U_e[i],
                ecolor=dp2_bar_color,
                linewidth=dp2_bar_line_width,
            )
            eb2[-1][0].set_linestyle(dp2_bar_line_style)
        elif i in ind:
            eb1 = axb.errorbar(
                Tall[i],
                Th_U[i],
                xerr=stats.norm.ppf(cr + (1 - cr) / 2.0) * s1[i],
                yerr=Th_U_e[i],
                ecolor=dp1_bar_color,
                linewidth=dp1_bar_line_width,
            )
            eb1[-1][0].set_linestyle(dp1_bar_line_style)

    # plot data point
    if opt_data_point:
        plot_data_point(axb, Tall, Th_U, ind, outd, outd_disc)


# ------------------------------------------------
# Figure setup
def makefigures(pd):
    # pd = plot_diagrams
    # ax1 = fig.add_subplot(221) # Conventional concordia
    # ax2 = fig.add_subplot(222) # Tera-Wasserburg concordia
    # ax3 = fig.add_subplot(223) # 238U/206Pb age
    # ax4 = fig.add_subplot(224) # Th/U - histogram

    if np.sum(pd) == 1:
        fig, ax = plt.subplots(2, 1, figsize=(6, 8))
        ax[1].axis("off")
    elif np.sum(pd) == 2:
        fig, ax = plt.subplots(1, 2, figsize=(12, 4))
    elif np.sum(pd) == 3:
        fig, ax = plt.subplots(3, 1, figsize=(6, 12))
    elif np.sum(pd) == 4:
        fig, ax = plt.subplots(2, 2, figsize=(12, 8))
    else:
        sys.exit("Error in plot_diagrams")

    mpl.rcParams["xtick.labelsize"] = legend_font_size
    mpl.rcParams["ytick.labelsize"] = legend_font_size
    ax = ax.ravel()

    if opt_correct_disequilibrium:
        figtitle = outfile + " (correct disequilibria"
        if opt_correct_common_Pb:
            figtitle = figtitle + " and correct common Pb"
    else:
        figtitle = outfile + " ("
        if opt_correct_common_Pb:
            figtitle = figtitle + "correct common Pb"
        else:
            figtitle = figtitle + "no correction"
    figtitle += ")"

    fig.canvas.set_window_title("%s" % figtitle)
    fig.suptitle("%s" % figtitle)

    return (fig, ax)


# plot KDE curves for all and accepted data
def plot_kde(ax_kde, rx, x, ii):
    # rx = range_hist_x
    # T = Tall
    # ii = ind
    ls = np.linspace(rx[0], rx[1], num=400)
    x = x[(x > rx[0]) & (x < rx[1])]
    if len(x) == 0:
        sys.exit("Please set appropriate axis age range in configuration file.")

    kde_all = stats.gaussian_kde(x)

    xi = x[ii]
    xi = xi[(xi > rx[0]) & (xi < rx[1])]
    kde = stats.gaussian_kde(xi)

    kde_multi_all = len(x)  # replace len(ls) 20190606
    kde_multi = len(xi)  # replace len(ls) 20190606

    ax_kde.plot(
        ls,
        kde_all(ls) * kde_multi_all,
        linestyle="--",
        color=kde_line_color,
        linewidth=kde_line_width,
    )
    ax_kde.plot(
        ls,
        kde(ls) * kde_multi,
        linestyle="-",
        color=kde_line_color,
        linewidth=kde_line_width,
    )


# plot histograms
def plot_hist(ax_hist, T, ii, od, oo):
    # T = Tall
    # ii = ind
    # od = outd_disc
    # oo = outd
    ax_hist.set_ylabel("Number of samples", fontsize=legend_font_size + 4)
    n, bins, rects = ax_hist.hist(
        (T[ii], T[od], T[oo]),
        hist_bin_num,
        histtype="barstacked",
        color=(hist_bin_color1, hist_bin_color2, hist_bin_color0),
        alpha=hist_bin_alpha,
        edgecolor="k",
        zorder=0,
        range=ax_hist.get_xlim(),
    )

    return (n, bins, rects)


# ################################################
if __name__ == "__main__":

    parser = OptionParser()
    parser.add_option(
        "-i",
        "--in",
        dest="inputfile",
        help="Name of input data file",
        metavar="FILE",
        action="store",
        type="string",
    )
    parser.add_option(
        "-c",
        "--cfg",
        dest="cfgfile",
        help="Name of configuration file",
        metavar="FILE",
        action="store",
        type="string",
    )
    parser.add_option(
        "-o",
        "--out",
        dest="outfile",
        help="Name of output file (when pdf driver is used)",
        metavar="FILE",
        action="store",
        type="string",
    )
    parser.add_option(
        "-d",
        "--driver",
        default="pdf",
        dest="driver",
        choices=["Qt5Agg", "qt5agg", "pdf", "TKAgg", "tkagg", "macosx"],
        help="Choose driver [Qt5Agg, TKAgg, maxosx, pdf (default)]",
    )
    parser.add_option(
        "-f",
        "--force-overwrite",
        help="Force overwrite the pre-existing pdf",
        default=False,
        action="store_true",
        dest="opt_force_overwrite",
    )
    (options, args) = parser.parse_args()

    mpl.use(options.driver)
    mpl.rcParams["pdf.fonttype"] = 42
    mpl.rcParams["ps.useafm"] = True
    mpl.rcParams["font.family"] = "Arial"
    import matplotlib.pyplot as plt

    # input filename is command-line argument or standard input
    if options.inputfile:
        infile = set_filename_input(options.inputfile)
    else:
        infile = set_filename_input()

    if options.cfgfile:
        conffile = set_filename_conf(options.cfgfile)
    else:
        conffile = set_filename_conf(infile)

    if options.outfile:
        outfile = set_filename_output(
            options.outfile, options.driver, options.opt_force_overwrite
        )
    else:
        outfile = set_filename_output(
            conffile, options.driver, options.opt_force_overwrite
        )

    # ################################################
    # Configuration
    config = ConfigParser(
        defaults={
            "infile_delimeter": "comma",
            "rows_of_header": "1",
            "colnum_207Pb_235U": "2",
            "colnum_207Pb_235U_error": "3",
            "colnum_206Pb_238U": "4",
            "colnum_206Pb_238U_error": "5",
            "colnum_207Pb_206Pb": "6",
            "colnum_207Pb_206Pb_error": "7",
            "colnum_207Pb_206Pb_inverse": False,
            "error_real": True,
            "input_error_sigma": "2",
            "opt_exclude_discordant_data": True,
            "discordance_percent_threshold": "10.0",
            "disc_type": "5",
            "opt_outlier": True,
            "outlier_alpha": "0.05",
            "exclude_data_points": "[]",
            "opt_Th_U": False,
            "opt_correct_disequilibrium": True,
            "f_Th_U": "0.2",
            "f_Pa_U": "3.5",
            "opt_correct_common_Pb": False,
            "Th_U_inverse": False,
            "Th_U_row_num": "[8]",
            "Th_U_error_num": "[]",
            "digits_number_output": "1",
            "legend_font_size": "8",
            "plot_diagrams": "[1, 1, 1, 1]",
            "opt_data_point": False,
            "opt_data_point_ee": True,
            "dp_ee_sigma": "2.0",
            "opt_2D_weighted_mean": True,
            "twm_ee_sigma": "2",
            "opt_concordia_age": True,
            "ca_sigma": "2",
            "ca_mswd": "1",
            "opt_concordia_intercept_age": False,
            "age_unit_name": "Ma",
            "graph_label_interval": "50",
            "graph_age_min": "5.0",
            "graph_age_max": "1500.0",
            "range_automatic_cc": True,
            "range_xy_cc": "[[0,0.05],[0.0,0.01]]",
            "range_automatic_twc": True,
            "range_xy_tw": "[[100, 1200], [0.00, 1.0]]",
            "range_automatic_oneD": True,
            "range_oneD_y": "[1, 110]",
            "oneD_age_type": "68",
            "oneD_yaxis_log": False,
            "oneD_plot_sort": True,
            "oneD_sigma": "2",
            "range_automatic_hist": True,
            "range_hist_x": "[0,1500]",
            "hist_age_type": "68",
            "Th_U_sigma": "2",
            "opt_kde": True,
            "dp0_marker_type": "s",
            "dp1_marker_type": "o",
            "dp2_marker_type": "^",
            "dp0_marker_size": "2",
            "dp1_marker_size": "2",
            "dp2_marker_size": "2",
            "dp0_marker_alpha": "1",
            "dp1_marker_alpha": "1",
            "dp2_marker_alpha": "1",
            "dp0_marker_face_color": "blue",
            "dp1_marker_face_color": "black",
            "dp2_marker_face_color": "red",
            "dp0_marker_edge_color": "blue",
            "dp1_marker_edge_color": "black",
            "dp2_marker_edge_color": "red",
            "dp0_marker_edge_width": "0.5",
            "dp1_marker_edge_width": "0.5",
            "dp2_marker_edge_width": "0.5",
            "dp0_ee_alpha": "1.0",
            "dp1_ee_alpha": "1.0",
            "dp2_ee_alpha": "1.0",
            "dp0_ee_face_color": "0.0",
            "dp1_ee_face_color": "0.0",
            "dp2_ee_face_color": "0.0",
            "dp0_ee_edge_line_style": "dotted",
            "dp1_ee_edge_line_style": "solid",
            "dp2_ee_edge_line_style": "dashed",
            "dp0_ee_edge_color": "blue",
            "dp1_ee_edge_color": "black",
            "dp2_ee_edge_color": "red",
            "dp0_ee_edge_width": "0.5",
            "dp1_ee_edge_width": "0.5",
            "dp2_ee_edge_width": "0.5",
            "twm_ee_alpha": "0.5",
            "twm_ee_face_color": "blue",
            "twm_ee_edge_color": "blue",
            "twm_ee_edge_width": "0.0",
            "ca_marker_type": "s",
            "ca_marker_size": "2",
            "ca_marker_face_color": "black",
            "ca_marker_edge_color": "black",
            "ca_marker_edge_width": "0.5",
            "concordia_ia_case_cc": "0",
            "concordia_ia_case_tw": "0",
            "ia_line_width": "1",
            "ia_line_color": "blue",
            "ia_sigma": "2",
            "ia_fill_color": "blue",
            "ia_alpha": "0.1",
            "oneD_wm_line_width": "2",
            "oneD_wm_line_color": "blue",
            "oned_band_fillcolor": "0.8",
            "oneD_band_alpha": "0.75",
            "dp0_bar_line_style": "dotted",
            "dp1_bar_line_style": "solid",
            "dp2_bar_line_style": "dashed",
            "dp0_bar_line_width": "0.5",
            "dp1_bar_line_width": "1.0",
            "dp2_bar_line_width": "0.5",
            "dp0_bar_color": "blue",
            "dp1_bar_color": "black",
            "dp2_bar_color": "red",
            "hist_bin_num": "20",
            "hist_bin_color0": "white",
            "hist_bin_color1": "blue",
            "hist_bin_color2": "0.5",
            "hist_bin_alpha": "0.5",
            "kde_line_color": "green",
            "kde_line_width": "1.5",
        }
    )

    config.read(conffile)

    c_delim = config.get("File", "infile_delimeter")  # 'comma' or 'tab'
    rows_of_header = loads(config.get("File", "rows_of_header"))
    c_7Pb5U = config.getint("File", "colnum_207Pb_235U")
    c_7Pb5U_e = config.getint("File", "colnum_207Pb_235U_error")
    c_6Pb8U = config.getint("File", "colnum_206Pb_238U")
    c_6Pb8U_e = config.getint("File", "colnum_206Pb_238U_error")
    c_7Pb6Pb = config.getint("File", "colnum_207Pb_206Pb")
    c_7Pb6Pb_e = config.getint("File", "colnum_207Pb_206Pb_error")
    c_7Pb6Pb_i = config.getboolean("File", "colnum_207Pb_206Pb_inverse")
    error_real = config.getboolean("File", "error_real")
    input_error_sigma = config.getfloat("File", "input_error_sigma")  # 1 or 2
    opt_exclude_disc = config.getboolean("File", "opt_exclude_discordant_data")
    disc_thres = config.getfloat("File", "discordance_percent_threshold")
    disc_type = config.getint("File", "disc_type")
    opt_outlier = config.getboolean("File", "opt_outlier")  # exclude outlier
    outlier_alpha = config.getfloat("File", "outlier_alpha")  # significant level
    exclude_data_points = loads(config.get("File", "exclude_data_points"))
    opt_Th_U = config.getboolean("File", "opt_Th_U")
    opt_correct_disequilibrium = config.getboolean("File", "opt_correct_disequilibrium")
    f_Th_U = config.getfloat("File", "f_Th_U")
    f_Pa_U = config.getfloat("File", "f_Pa_U")
    opt_correct_common_Pb = config.getboolean("File", "opt_correct_common_Pb")
    Th_U_inverse = config.getboolean("File", "Th_U_inverse")
    Th_U_row_num = loads(config.get("File", "Th_U_row_num"))
    Th_U_error_num = loads(config.get("File", "Th_U_error_num"))
    dig_num_output = config.getint("Graph", "digits_number_output")  # 2
    plot_diagrams = loads(config.get("Graph", "plot_diagrams"))  # [1, 1, 1, 1]
    graph_age_min = config.getfloat("Graph", "graph_age_min")
    graph_age_max = config.getfloat("Graph", "graph_age_max")
    graph_label_interval = config.getfloat("Graph", "graph_label_interval")
    age_unit_name = config.get("Graph", "age_unit_name")  # = 'Ma'
    legend_font_size = config.getint("Graph", "legend_font_size")  # = 10
    range_automatic_cc = config.getboolean("Graph", "range_automatic_cc")
    range_xy_cc = loads(config.get("Graph", "range_xy_cc"))  # [[0,6],[0,0.35]]
    range_automatic_twc = config.getboolean("Graph", "range_automatic_twc")
    range_xy_tw = loads(config.get("Graph", "range_xy_tw"))  # [[xmin,xmax],[ymin,ymax]]
    opt_data_point = config.getboolean("Graph", "opt_data_point")
    dp0_marker_type = config.get("Graph", "dp0_marker_type")  # = 'o'
    dp0_marker_size = config.getfloat("Graph", "dp0_marker_size")  # = 7
    dp0_marker_alpha = config.getfloat("Graph", "dp0_marker_alpha")  # = 1.0
    dp0_marker_fc = config.get("Graph", "dp0_marker_face_color")  # white
    dp0_marker_ec = config.get("Graph", "dp0_marker_edge_color")  # black
    dp0_marker_ew = config.getfloat("Graph", "dp0_marker_edge_width")  # 0.5
    dp1_marker_type = config.get("Graph", "dp1_marker_type")  # = 'o'
    dp1_marker_size = config.getfloat("Graph", "dp1_marker_size")  # = 7
    dp1_marker_alpha = config.getfloat("Graph", "dp1_marker_alpha")  # = 1.0
    dp1_marker_fc = config.get("Graph", "dp1_marker_face_color")  # black
    dp1_marker_ec = config.get("Graph", "dp1_marker_edge_color")  # white
    dp1_marker_ew = config.getfloat("Graph", "dp1_marker_edge_width")  # 0.5
    dp2_marker_type = config.get("Graph", "dp2_marker_type")
    dp2_marker_size = config.getfloat("Graph", "dp2_marker_size")
    dp2_marker_alpha = config.getfloat("Graph", "dp2_marker_alpha")
    dp2_marker_fc = config.get("Graph", "dp2_marker_face_color")
    dp2_marker_ec = config.get("Graph", "dp2_marker_edge_color")
    dp2_marker_ew = config.getfloat("Graph", "dp2_marker_edge_width")
    opt_data_point_ee = config.getboolean("Graph", "opt_data_point_ee")
    dp_ee_sigma = config.getfloat("Graph", "dp_ee_sigma")  # 2
    dp0_ee_alpha = config.getfloat("Graph", "dp0_ee_alpha")  # 1.0
    dp0_ee_fc = config.get("Graph", "dp0_ee_face_color")  # none
    dp0_ee_ls = config.get("Graph", "dp0_ee_edge_line_style")  # '-'
    dp0_ee_ec = config.get("Graph", "dp0_ee_edge_color")  # 0.5
    dp0_ee_ew = config.get("Graph", "dp0_ee_edge_width")  # 0.5
    dp1_ee_alpha = config.getfloat("Graph", "dp1_ee_alpha")  # 0.2
    dp1_ee_fc = config.get("Graph", "dp1_ee_face_color")  # 0.5
    dp1_ee_ls = config.get("Graph", "dp1_ee_edge_line_style")  # '-'
    dp1_ee_ec = config.get("Graph", "dp1_ee_edge_color")  # white
    dp1_ee_ew = config.get("Graph", "dp1_ee_edge_width")  # 0.5
    dp2_ee_alpha = config.getfloat("Graph", "dp2_ee_alpha")
    dp2_ee_fc = config.get("Graph", "dp2_ee_face_color")
    dp2_ee_ls = config.get("Graph", "dp2_ee_edge_line_style")  # ':'
    dp2_ee_ec = config.get("Graph", "dp2_ee_edge_color")
    dp2_ee_ew = config.get("Graph", "dp2_ee_edge_width")
    opt_2D_wm = config.getboolean("Graph", "opt_2D_weighted_mean")
    twm_ee_sigma = config.getfloat("Graph", "twm_ee_sigma")  # 2
    twm_ee_fc = config.get("Graph", "twm_ee_face_color")  # green
    twm_ee_ec = config.get("Graph", "twm_ee_edge_color")  # none
    twm_ee_ew = config.get("Graph", "twm_ee_edge_width")  # 0.5
    twm_ee_alpha = config.getfloat("Graph", "twm_ee_alpha")  # 0.7
    opt_concordia_age = config.getboolean("Graph", "opt_concordia_age")
    concordia_ia_case_cc = config.getint("Graph", "concordia_ia_case_cc")
    concordia_ia_case_tw = config.getint("Graph", "concordia_ia_case_tw")
    ca_sigma = config.getfloat("Graph", "ca_sigma")  # 2
    ca_marker_type = config.get("Graph", "ca_marker_type")  # s
    ca_marker_size = config.getfloat("Graph", "ca_marker_size")  # 8
    ca_marker_fc = config.get("Graph", "ca_marker_face_color")  # magenta
    ca_marker_ec = config.get("Graph", "ca_marker_edge_color")  # black
    ca_marker_ew = config.getfloat("Graph", "ca_marker_edge_width")  # 1.0
    ca_mswd = config.getint("Graph", "ca_mswd")  # 0, 1, or 2
    opt_concordia_ia = config.getboolean("Graph", "opt_concordia_intercept_age")
    ia_line_width = config.getfloat("Graph", "ia_line_width")  # 1
    ia_line_color = config.get("Graph", "ia_line_color")  # blue
    ia_sigma = config.getfloat("Graph", "ia_sigma")  # 2
    ia_fill_color = config.get("Graph", "ia_fill_color")  # blue
    ia_alpha = config.getfloat("Graph", "ia_alpha")  # 0.1
    range_automatic_oneD = config.getboolean("Graph", "range_automatic_oneD")
    range_oneD_y = loads(config.get("Graph", "range_oneD_y"))  # [70, 100]
    oneD_age_type = config.getint("Graph", "oneD_age_type")  # 68
    oneD_plot_sort = config.getboolean("Graph", "oneD_plot_sort")
    oneD_sigma = config.getfloat("Graph", "oneD_sigma")  # 2
    oneD_wm_line_width = config.getfloat("Graph", "oneD_wm_line_width")  # 2
    oneD_wm_line_color = config.get("Graph", "oneD_wm_line_color")  # blue
    oneD_band_fc = config.get("Graph", "oneD_band_fillcolor")  # 0.8
    oneD_band_alpha = config.getfloat("Graph", "oneD_band_alpha")  # 0.5
    oneD_yaxis_log = config.getboolean("Graph", "oneD_yaxis_log")  # 0
    dp0_bar_line_style = config.get("Graph", "dp0_bar_line_style")  # solid
    dp0_bar_line_width = config.getfloat("Graph", "dp0_bar_line_width")  # 1
    dp0_bar_color = config.get("Graph", "dp0_bar_color")  # black
    dp1_bar_line_style = config.get("Graph", "dp1_bar_line_style")  # solid
    dp1_bar_line_width = config.getfloat("Graph", "dp1_bar_line_width")  # 1
    dp1_bar_color = config.get("Graph", "dp1_bar_color")  # black
    dp2_bar_line_style = config.get("Graph", "dp2_bar_line_style")  # dashed
    dp2_bar_line_width = config.getfloat("Graph", "dp2_bar_line_width")  # 1
    dp2_bar_color = config.get("Graph", "dp2_bar_color")  # black
    range_automatic_hist = config.getboolean("Graph", "range_automatic_hist")
    range_hist_x = loads(config.get("Graph", "range_hist_x"))  # [70, 100]
    hist_bin_num = config.getint("Graph", "hist_bin_num")  # 20
    hist_age_type = config.getint("Graph", "hist_age_type")  # 68
    Th_U_sigma = config.getfloat("Graph", "Th_U_sigma")  # 2
    hist_bin_color0 = config.get("Graph", "hist_bin_color0")  # white
    hist_bin_color1 = config.get("Graph", "hist_bin_color1")  # blue
    hist_bin_color2 = config.get("Graph", "hist_bin_color2")  # 0.5
    hist_bin_alpha = config.getfloat("Graph", "hist_bin_alpha")  # 0.75
    opt_kde = config.getboolean("Graph", "opt_kde")  # 1
    kde_line_color = config.get("Graph", "kde_line_color")  # red
    kde_line_width = config.get("Graph", "kde_line_width")  # 2

    # cumulative probability density
    # 1 sigma (68.3), 2 sigma (95.4%), and 3 sigma (99.7%)
    dp_ee_cr = 2.0 * stats.norm.cdf(dp_ee_sigma) - 1.0
    twm_ee_cr = 2.0 * stats.norm.cdf(twm_ee_sigma) - 1.0
    ca_cr = 2.0 * stats.norm.cdf(ca_sigma) - 1.0
    ia_cr = 2.0 * stats.norm.cdf(ia_sigma) - 1.0
    oneD_cr = 2.0 * stats.norm.cdf(oneD_sigma) - 1.0
    Th_U_cr = 2.0 * stats.norm.cdf(Th_U_sigma) - 1.0

    # Define age unit
    if age_unit_name == "ka":
        age_unit = 10 ** 3
        time = time_ka
    else:
        age_unit = 10 ** 6
        time = time_ma

    # ################################################
    # Data formatting

    # format of output digits number
    dignum = "." + str(dig_num_output) + "f"
    dignum2 = "." + str(2) + "f"

    # ------------------------------------------------
    # Set automatic ranges for plotting

    # Conventional concordia diagrams
    if range_automatic_cc:
        range_XY = [
            [
                np.exp(l235U * graph_age_min * age_unit) - 1,
                np.exp(l235U * graph_age_max * age_unit) - 1,
            ],
            [
                np.exp(l238U * graph_age_min * age_unit) - 1,
                np.exp(l238U * graph_age_max * age_unit) - 1,
            ],
        ]
    else:
        range_XY = range_xy_cc

    rX = [range_XY[0][0], range_XY[0][1]]
    rY = [range_XY[1][0], range_XY[1][1]]
    tX_min, tX_max, tY_min, tY_max = TimeRangeConv(rX, rY)
    timeX = [t for t in time if t >= tX_min and t <= tX_max]
    timeY = [t for t in time if t >= tY_min and t <= tY_max]
    if len(timeX) < len(timeY):
        timeXY = timeX + [tX_max * 1.1]
        timeXY.insert(0, tX_min)
    else:
        timeXY = timeY + [tY_max * 1.1]
        timeXY.insert(0, tY_min)

    Xconv, Yconv = ConcLineConv(np.array(timeXY))
    line_cc = [[Xconv[i], Yconv[i]] for i, j in enumerate(Xconv)]

    # ------------------------------------------------
    # Tera-Wasserburg concordia diagrams
    # x = 1/Y
    # y = 1/137.82*X/Y
    if range_automatic_twc:
        range_xy = [
            [1.0 / range_XY[1][1], 1.0 / range_XY[1][0]],
            [
                1.0 / U85r * range_XY[0][0] / range_XY[1][0],
                1.0 / U85r * range_XY[0][1] / range_XY[1][1],
            ],
        ]
    else:
        range_xy = range_xy_tw

    rx = [range_xy[0][0], range_xy[0][1]]
    ry = [range_xy[1][0], range_xy[1][1]]
    tx_min, tx_max = TimeRangeTW(rx)
    timexy = [float(t) for t in time if t >= tx_min and t <= tx_max]
    timexy += [tx_max * 1.1]
    timexy.insert(0, tx_min)
    Xtw, Ytw = ConcLineTW(np.array(timexy))
    line_tw = [[Xtw[i], Ytw[i]] for i, j in enumerate(Xtw)]

    # 1D bar plot
    # x = number of samples
    # y = 206Pb/238U, 207Pb/235U, or 207Pb/206Pb age
    if range_automatic_oneD:
        range_oneD_y = [graph_age_min, graph_age_max]

    # histogram
    # y1 = number of samples
    # y2 = Th/U (optional)
    # x = 206Pb/238U, 207Pb/235U, or 207Pb/206Pb age
    if range_automatic_hist:
        range_hist_x = [graph_age_min, graph_age_max]

    # ================================================
    # Data formatting

    dt_name_column = [
        "r75",
        "r75e",
        "r68",
        "r68e",
        "r76",
        "r76e",
    ]

    column_num_isotopic_ratio = [
        c_7Pb5U,
        c_7Pb5U_e,
        c_6Pb8U,
        c_6Pb8U_e,
        c_7Pb6Pb,
        c_7Pb6Pb_e,
    ]

    if c_delim == "comma":
        delim = ","
    elif c_delim == "tab":
        delim = "\t"
    else:
        delim = " "

    d = np.loadtxt(
        infile,
        delimiter=delim,
        usecols=column_num_isotopic_ratio,
        skiprows=rows_of_header,
    )

    d = pd.DataFrame(d, columns=dt_name_column)

    # if 207Pb/206Pb is given by 206Pb/207Pb
    if c_7Pb6Pb_i:
        d["r76"] = 1.0 / d["r76"]

    # if error is shown in percentage
    if not error_real:
        d["r75e"] = d["r75"] * d["r75e"] / 100.0
        d["r68e"] = d["r68"] * d["r68e"] / 100.0
        d["r76e"] = d["r76"] * d["r76e"] / 100.0
    # if error range is given by other than 1 sigma
    if input_error_sigma != 1.0:
        d["r75e"] = d["r75e"] / input_error_sigma
        d["r68e"] = d["r68e"] / input_error_sigma
        d["r76e"] = d["r76e"] / input_error_sigma

    # Conventional concordia diagrams
    # (X, Y) = (207Pb/235U, 206Pb/238U)
    # Conventional concordia diagram
    X = d["r75"]
    Y = d["r68"]

    sigma_X = d["r75e"]
    sigma_Y = d["r68e"]
    SX = sigma_X / X
    SY = sigma_Y / Y

    # Tera-Wasserburg concordia diagrams
    # (x, y) = (238U/206Pb, 207Pb/206Pb)
    x = 1 / Y
    y = d["r76"]

    sigma_x = SY * x
    sigma_y = d["r76e"]
    # Sx = sigma_x/x
    Sx = SY
    Sy = sigma_y / y

    # error correlation
    rho_XY = (SX ** 2 + SY ** 2 - Sy ** 2) / (2.0 * SX * SY)
    # rho_xy = (SY**2-SX**2*rho_XY)/Sy # Equation in p. 27 of Ludwig2012
    rho_xy = (SY ** 2 - SX * SY * rho_XY) / (Sx * Sy)

    # covariance
    cov_XY = rho_XY * sigma_X * sigma_Y
    cov_xy = rho_xy * sigma_x * sigma_y

    if opt_Th_U:
        Th_U = np.loadtxt(
            infile, delimiter=delim, usecols=Th_U_row_num, skiprows=rows_of_header
        )
        if Th_U_inverse:
            Th_U = 1.0 / Th_U
        if Th_U_error_num:
            Th_U_e = np.loadtxt(
                infile, delimiter=delim, usecols=Th_U_error_num, skiprows=rows_of_header
            )
        else:
            Th_U_e = Th_U * 0.0

    # ------------------------------------------------
    # Common Pb correction (207Pb method)
    # (Williams, 1998)
    #
    if opt_correct_common_Pb:
        # Correction of common Pb (207Pb-corrected)
        # xcorr, ycorr, f206, T68corr = Pb76c(x, y)
        d["86cor"], d["76cor"], d["f206p"], d["t68cor"] = func_Pb76c(x, y)
        x = d["86cor"]
        y = d["76cor"]
        Y = 1 / x
        Sy = sigma_y / y
        Sx = sigma_x / x
        SY = sigma_Y / Y
        # rho_XY = (SX ** 2 + SY ** 2 - Sy ** 2) / (2.0 * SX * SY)
        # rho_xy = (SY ** 2 - SX * SY * rho_XY) / (Sx * Sy)
        # cov_XY = rho_XY * sigma_X * sigma_Y
        # cov_xy = rho_xy * sigma_x * sigma_y

    # ------------------------------------------------
    # Ages
    d["t75"] = 1 / l235U * np.log(X + 1)
    d["t75e"] = np.empty(len(y))
    d["t68"] = 1 / l238U * np.log(Y + 1)
    d["t68e"] = np.empty(len(y))
    d["t76"] = np.empty(len(y))
    d["t76e_plus"] = np.empty(len(y))  # 1sigma error
    d["t76e_minus"] = np.empty(len(y))  # 1sigma error
    d["t76_min"] = np.empty(len(y))  # 1sigma error
    d["t76_max"] = np.empty(len(y))  # 1sigma error
    (d["t76"], d["t76e_plus"], d["t76e_minus"]) = calc_t76(age_unit, y, sigma_y, ca_cr)
    if opt_correct_disequilibrium:
        d["t68"] = SI_Tdiseq(d["r75"], d["r68"], rtype=68)
        d["t75"] = SI_Tdiseq(d["r75"], d["r68"], rtype=75)

    # ------------------------------------------------
    # Initialize list
    ind = np.array(range(len(x)))  # inliers
    outd_disc = []  # discordants
    outd = []  # outliers

    # ------------------------------------------------
    # Data points of additional exclusion
    excluded_points = [i for i in exclude_data_points if i is not None]
    if excluded_points:
        outd = np.unique(excluded_points)
        ind = np.delete(ind, outd)

    # ------------------------------------------------
    # Exclusion of discordant data for calculation
    if opt_exclude_disc:
        if disc_type == 5:
            # 206Pb/238U--207Pb/235U
            ind, outd_disc = discordant_judge(
                X, Y, sigma_X, sigma_Y, cov_XY, dp_ee_cr, ind, outd_disc, line_cc
            )
            # 207Pb/206Pb--238U/206Pb
            ind, outd_disc = discordant_judge(
                x, y, sigma_x, sigma_y, cov_xy, dp_ee_cr, ind, outd_disc, line_tw
            )
        else:
            disc_percent = discordance(
                d["t75"],
                d["t75e"],
                d["t68"],
                d["t68e"],
                d["t76"],
                d["t76_min"],
                d["t76_max"],
                input_error_sigma,
                method=disc_type,
            )
            # outd_disc = np.where((np.abs(disc_percent) >= disc_thres) | (disc_percent < 0))
            if np.max(disc_percent) >= disc_thres:
                outd_disc = np.where(np.abs(disc_percent) >= disc_thres)
    else:
        outd_disc = []
        ind = np.delete(ind, outd_disc)

    # ################################################
    # List of the configurations

    # input file
    print("\n")
    print("============================================================")
    print(("Data filename is %s") % infile)
    print(("Configuration filename is %s") % conffile)
    print(("Output filename is %s") % outfile)
    print("# Input data (first 5 lines)")
    print("------------------------------------------------------------")
    print("207Pb/235U  1s     206Pb/238U  1s     207Pb/206Pb  1s")
    print(
        "column[%d]   column[%d] column[%d]   column[%d] column[%d]    column[%d]"
        % (
            column_num_isotopic_ratio[0],
            column_num_isotopic_ratio[1],
            column_num_isotopic_ratio[2],
            column_num_isotopic_ratio[3],
            column_num_isotopic_ratio[4],
            column_num_isotopic_ratio[5],
        )
    )

    if len(X) > 5:
        Nc = 5
    else:
        Nc = len(X)

    for i in range(0, Nc):
        print(
            "%.5f     %.5f   %.5f    %.5f    %.5f    %.5f"
            % (X[i], sigma_X[i], Y[i], sigma_Y[i], y[i], sigma_y[i])
        )

    print("------------------------------------------------------------")
    # Range of data
    print("Range of data")
    print("207Pb/235U: %.5f--%.5f" % (d["r75"].min(), d["r75"].max()))
    print("206Pb/238U: %.5f--%.5f" % (d["r68"].min(), d["r68"].max()))
    print("207Pb/206Pb: %.5f--%.5f" % (d["r76"].min(), d["r76"].max()))
    print("238U/206Pb: %.5f--%.5f" % (1 / d["r68"].max(), 1 / d["r68"].min()))

    print("------------------------------------------------------------")
    # discordance
    if opt_exclude_disc:
        if disc_type == 5:
            print(
                "Discordant data means that the error ellipses do not intersect the concordia line."
            )
        else:
            print_discordant_data(disc_type, outd_disc, disc_percent, input_error_sigma)
    else:
        print("Discordant data are not excluded from calculation")

    # ------------------------------------------------
    # exit if accepted data < 2
    if len(ind) < 2:
        sys.exit("Error: accepted data points are less than 2.")

    # Check correlation coefficient
    if rho_XY[ind].max() > 1:
        print("Correlation coefficient: rho_XY")
        [print("%d %.2f" % (i, rho_XY[i])) for i in range(len(rho_XY)) if rho_XY[i] > 1]
        sys.exit("rho_XY is more than 1")
    elif rho_xy[ind].max() > 1:
        print("Correlation coefficient: rho_xy")
        [print("%d %.2f" % (i, rho_xy[i])) for i in range(len(rho_xy)) if rho_xy[i] > 1]
        sys.exit("rho_xy is more than 1")

    # ------------------------------------------------
    # corrections for common Pb and initial disequilibria
    if opt_correct_common_Pb:
        if opt_correct_disequilibrium:
            d["t68"] = d["t68cor"]
        else:
            d["t68"] = (1 - d["f206p"] / 100) * d["t68cor"]

    # ################################################
    # print ages
    print("------------------------------------------------------------")
    if opt_correct_disequilibrium:
        print("U-Pb ages^ (%s) [%d sigma]" % (age_unit_name, ca_sigma))
    else:
        print("U-Pb ages (%s) [%d sigma]" % (age_unit_name, ca_sigma))

    if opt_correct_common_Pb:
        print("#\tf206%\t6/8*\t+-s\t7/5\t+-s\t7/6*\t+s\t-s\t")
    else:
        print("#\t6/8\t+-s\t7/5\t+-s\t7/6\t+s\t-s\t")

    for i in range(len(y)):
        d["t68e"].iloc[i] = d["t68"].iloc[i] * SY[i]
        d["t75e"].iloc[i] = d["t75"].iloc[i] * SX[i]

        if opt_correct_common_Pb:
            print(
                "%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s"
                % (
                    i + 1,
                    format(d["f206p"].iloc[i], dignum),
                    format(d["t68"].iloc[i] / age_unit, dignum),
                    format(d["t68e"].iloc[i] / age_unit, dignum),
                    format(d["t75"].iloc[i] / age_unit, dignum),
                    format(d["t75e"].iloc[i] / age_unit, dignum),
                    format(d["t76"].iloc[i] / age_unit, dignum),
                    format(d["t76e_plus"].iloc[i] / age_unit, dignum),
                    format(d["t76e_minus"].iloc[i] / age_unit, dignum),
                )
            )
        else:
            print(
                "%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s"
                % (
                    i + 1,
                    format(d["t68"].iloc[i] / age_unit, dignum),
                    format(d["t68e"].iloc[i] / age_unit, dignum),
                    format(d["t75"].iloc[i] / age_unit, dignum),
                    format(d["t75e"].iloc[i] / age_unit, dignum),
                    format(d["t76"].iloc[i] / age_unit, dignum),
                    format(d["t76e_plus"].iloc[i] / age_unit, dignum),
                    format(d["t76e_minus"].iloc[i] / age_unit, dignum),
                )
            )

    if opt_correct_common_Pb:
        print("* Correction for common Pb by 207Pb method")
    if opt_correct_disequilibrium:
        print("^ Correction for initial disequilibria.")

    # ################################################
    print("------------------------------------------------------------")
    # Excluding outliers
    if opt_outlier:
        Tall, s1, label_selected = select_age_type(oneD_age_type)

        # check for normality before GESD test
        w, p = stats.shapiro(list(Tall))
        if p < outlier_alpha:
            print("Normality test: passed (p = %.2f)" % (p))
        else:
            print("Normality test: failed (p = %.2f)" % (p))

        # generalize ESD test
        ii, oo = GESDtest(Tall, s1, ind, outlier_alpha)

        if len(ii) > 0:
            ind = ii
            print("Concordants (accepted) are ", end=" ")
            print(ind + 1)
        if len(oo) > 0:
            outd = np.unique(np.append(outd, oo))
            print("Concordants (excluded) are ", end=" ")
            print(outd + 1)
        if len(outd_disc) > 0:
            print("Discordants are ", end=" ")
            print(outd_disc + 1)
    else:
        print("Concordants (accepted) are ", end=" ")
        print(ind)
        print("Concordants (excluded) are ", end=" ")
        print(outd)
        print("Discordants are ", end=" ")
        print(outd_disc)

    # excluded data points
    if excluded_points:
        # print('Manually excluded data points are'),  # python2
        print("Manually excluded data points are", end=" ")  # python3
        print(excluded_points)

    # ------------------------------------------------
    # Number of data points

    # total number
    N = len(X)

    # ################################################
    # plotting
    fig, ax = makefigures(plot_diagrams)

    # ------------------------------------------------
    # A: Conventional concordia plot

    if plot_diagrams[0] == 1:
        axn = 0
        axn_title = "a"

        print("------------------------------------------------------------")
        print(("%s: Conventional concordia diagram") % axn_title)

        ax[axn].set_title(axn_title, loc="left", fontsize=legend_font_size + 6)
        ax[axn].set_xlabel("$^{207}$Pb* / $^{235}$U", fontsize=legend_font_size + 4)
        ax[axn].set_ylabel("$^{206}$Pb* / $^{238}$U", fontsize=legend_font_size + 4)
        ax[axn].set_xlim(rX)
        ax[axn].set_ylim(rY)

        PlotConcConv(
            ax[axn],
            Xconv,
            Yconv,
            timeXY,
            age_unit,
            math.ceil(graph_label_interval),
            legend_font_size,
        )

        # Legend
        legend_pos_x, legend_pos_y = calc_legend_pos(range_XY)
        legend_pos = 0

        # draw error ellipses
        if opt_data_point_ee:
            print("    Error ellipses are %d%% for data points" % (dp_ee_cr * 100))

            if disc_type == 5:
                line_cc = [[Xconv[i], Yconv[i]] for i, j in enumerate(Xconv)]
                ind_cc, outd_disc_cc = plot_data_point_error_ellipse_disc(
                    ax,
                    axn,
                    X,
                    Y,
                    sigma_X,
                    sigma_Y,
                    cov_XY,
                    dp_ee_cr,
                    ind,
                    outd,
                    outd_disc,
                    line_cc,
                )
                outd_disc = np.unique(np.append(outd_disc, outd_disc_cc))
                if len(ind_cc) > 0:
                    ind = ind_cc
                if len(outd_disc_cc) > 0:
                    outd_disc = outd_disc_cc
            else:
                plot_data_point_error_ellipse(
                    ax,
                    axn,
                    X,
                    Y,
                    sigma_X,
                    sigma_Y,
                    cov_XY,
                    dp_ee_cr,
                    ind,
                    outd,
                    outd_disc,
                )

        # plot data point
        if opt_data_point:
            plot_data_point(ax[axn], X, Y, ind, outd, outd_disc)

        # # weighted mean
        if opt_2D_wm:
            print("    Error ellipse is %d%% for 2D weighted mean" % (twm_ee_cr * 100))
            legend_pos += 1
            plot_2D_wm(
                ax,
                axn,
                X[ind],
                Y[ind],
                sigma_X[ind],
                sigma_Y[ind],
                rho_XY[ind],
                twm_ee_cr,
                legend_pos_x[legend_pos],
                legend_pos_y[legend_pos],
            )

        # Concordia age
        if opt_concordia_age or opt_concordia_ia:
            (
                T_lsq,
                S_lsq,
                X_lsq,
                Y_lsq,
                MSWDconc,
                MSWDeq,
                MSWDcomb,
                Pconc,
                Peq,
                Pcomb,
            ) = concordia_age(
                "conv", X[ind], Y[ind], sigma_X[ind], sigma_Y[ind], rho_XY[ind], ca_cr
            )

            if ca_mswd == 0:
                MSWD = MSWDconc
                Pvalue = Pconc
            elif ca_mswd == 1:
                MSWD = MSWDeq
                Pvalue = Peq
            else:
                MSWD = MSWDcomb
                Pvalue = Pcomb

        if opt_concordia_age:
            legend_pos += 1
            plot_concordia_age(
                ax,
                axn,
                T_lsq,
                S_lsq,
                X_lsq,
                Y_lsq,
                MSWD,
                ca_cr,
                legend_pos_x[legend_pos],
                legend_pos_y[legend_pos],
            )
            legend_pos += 1
            plot_concordia_age_MSWD(
                ax,
                axn,
                MSWD,
                ca_mswd,
                Pvalue,
                legend_pos_x[legend_pos],
                legend_pos_y[legend_pos],
            )

            print(
                u"    Concordia age = %s ± %s [%d%% conf.] / ± %s [t√MSWD] %s"
                % (
                    format(T_lsq / age_unit, dignum),
                    format(S_lsq / age_unit, dignum),
                    (ca_cr * 100),
                    format(S_lsq / age_unit * np.sqrt(MSWDcomb), dignum),
                    age_unit_name,
                )
            )
            print(
                "    MSWD concordance / equivalence / combined = %.2f / %.2f / %.2f"
                % (MSWDconc, MSWDeq, MSWDcomb)
            )
            print(
                "    P-value concordance / equivalence / combined = %.2f / %.2f / %.2f"
                % (Pconc, Peq, Pcomb)
            )

        # plot intercept line and band
        if opt_concordia_ia:
            if (concordia_ia_case_cc == 0) or (concordia_ia_case_cc == 2):
                legend_pos += 1
                case = 0
                Tsi, Tmax, Tmin = plot_concordia_intercept_age(
                    ax,
                    axn,
                    "conv",
                    X[ind],
                    Y[ind],
                    sigma_X[ind],
                    sigma_Y[ind],
                    ia_cr,
                    rho_XY[ind],
                    range_XY,
                    T_lsq,
                    case,
                    legend_pos_x[legend_pos],
                    legend_pos_y[legend_pos],
                )
                print(
                    "    Intercept age = %s +%s %s %s [%d%% conf.]"
                    % (
                        format(Tsi / age_unit, dignum),
                        format((Tmax - Tsi) / age_unit, dignum),
                        format((Tmin - Tsi) / age_unit, dignum),
                        age_unit_name,
                        ia_cr * 100,
                    )
                )

            if (concordia_ia_case_cc == 1) or (concordia_ia_case_cc == 2):
                legend_pos += 1
                case = 1
                Tsi, Tmax, Tmin = plot_concordia_intercept_age(
                    ax,
                    axn,
                    "conv",
                    X[ind],
                    Y[ind],
                    sigma_X[ind],
                    sigma_Y[ind],
                    ia_cr,
                    rho_XY[ind],
                    range_XY,
                    T_lsq,
                    case,
                    legend_pos_x[legend_pos],
                    legend_pos_y[legend_pos],
                )
                print(
                    "    Intercept age = %s +%s %s %s [%d%% conf.]"
                    % (
                        format(Tsi / age_unit, dignum),
                        format((Tmax - Tsi) / age_unit, dignum),
                        format((Tmin - Tsi) / age_unit, dignum),
                        age_unit_name,
                        ia_cr * 100,
                    )
                )

        # Sample number
        legend_data_number(
            ax, axn, legend_pos_x[legend_pos], legend_pos_y[legend_pos], ind
        )

    # ------------------------------------------------
    # B: Tera-Wasserburg concordia plot

    if plot_diagrams[1] == 1:
        if plot_diagrams[0] == 1:
            axn = 1
            axn_title = "b"
        else:
            axn = 0
            axn_title = "a"

        print(("%s: Tera-Wasserburg concordia diagram") % axn_title)

        ax[axn].set_title(axn_title, loc="left", fontsize=legend_font_size + 6)
        ax[axn].set_xlabel("$^{238}$U / $^{206}$Pb*", fontsize=legend_font_size + 4)
        ax[axn].set_ylabel("$^{207}$Pb* / $^{206}$Pb*", fontsize=legend_font_size + 4)
        ax[axn].set_xlim(rx)
        ax[axn].set_ylim(ry)

        PlotConcTW(
            ax[axn],
            Xtw,
            Ytw,
            timexy,
            age_unit,
            math.ceil(graph_label_interval),
            legend_font_size,
        )

        # Legend data number
        legend_pos_x, legend_pos_y = calc_legend_pos(range_xy)
        legend_pos = 0

        # plot data point
        if opt_data_point:
            plot_data_point(ax[axn], x, y, ind, outd, outd_disc)

        # draw error ellipses
        if opt_data_point_ee:
            if disc_type == 5:
                ind_tw, outd_disc_tw = plot_data_point_error_ellipse_disc(
                    ax,
                    axn,
                    x,
                    y,
                    sigma_x,
                    sigma_y,
                    cov_xy,
                    dp_ee_cr,
                    ind,
                    outd,
                    outd_disc,
                    line_tw,
                )
                outd_disc = np.unique(np.append(outd_disc, outd_disc_tw))
                if len(ind_tw) > 0:
                    ind = ind_tw
                if len(outd_disc_tw) > 0:
                    outd_disc = outd_disc_tw
            else:
                plot_data_point_error_ellipse(
                    ax,
                    axn,
                    x,
                    y,
                    sigma_x,
                    sigma_y,
                    cov_xy,
                    dp_ee_cr,
                    ind,
                    outd,
                    outd_disc,
                )

        # weighted mean
        if opt_2D_wm:
            legend_pos += 1
            plot_2D_wm(
                ax,
                axn,
                x[ind],
                y[ind],
                sigma_x[ind],
                sigma_y[ind],
                rho_xy[ind],
                twm_ee_cr,
                legend_pos_x[legend_pos],
                legend_pos_y[legend_pos],
            )

        # Concordia age
        if opt_concordia_age or opt_concordia_ia:
            (
                t_lsq,
                s_lsq,
                x_lsq,
                y_lsq,
                mswd_conc,
                mswd_eq,
                mswd_comb,
                p_conc,
                p_eq,
                p_comb,
            ) = concordia_age(
                "tw", x[ind], y[ind], sigma_x[ind], sigma_y[ind], rho_xy[ind], ca_cr
            )

            if ca_mswd == 0:
                mswd = mswd_conc
                pvalue = p_conc
            elif ca_mswd == 1:
                mswd = mswd_eq
                pvalue = p_eq
            else:
                mswd = mswd_comb
                pvalue = p_comb

        if opt_concordia_age:
            legend_pos += 1
            plot_concordia_age(
                ax,
                axn,
                t_lsq,
                s_lsq,
                x_lsq,
                y_lsq,
                mswd,
                ca_cr,
                legend_pos_x[legend_pos],
                legend_pos_y[legend_pos],
            )
            legend_pos += 1
            plot_concordia_age_MSWD(
                ax,
                axn,
                mswd,
                ca_mswd,
                pvalue,
                legend_pos_x[legend_pos],
                legend_pos_y[legend_pos],
            )

            print(
                u"    Concordia age = %s ± %s [%d%% conf.] / ± %s [t√MSWD] %s"
                % (
                    format(t_lsq / age_unit, dignum),
                    format(s_lsq / age_unit, dignum),
                    (ca_cr * 100),
                    format(s_lsq / age_unit * np.sqrt(mswd_comb), dignum),
                    age_unit_name,
                )
            )

            print(
                "    MSWD concordance / equivalence / combined = %.2f / %.2f / %.2f"
                % (mswd_conc, mswd_eq, mswd_comb)
            )
            print(
                "    P-value concordance / equivalence / combined = %.2f / %.2f / %.2f"
                % (p_conc, p_eq, p_comb)
            )

        # plot intercept line and band
        if opt_concordia_ia:
            if (concordia_ia_case_tw == 0) or (concordia_ia_case_tw == 2):
                legend_pos += 1
                case = 0
                Tsi, Tmax, Tmin = plot_concordia_intercept_age(
                    ax,
                    axn,
                    "tw",
                    x[ind],
                    y[ind],
                    sigma_x[ind],
                    sigma_y[ind],
                    ia_cr,
                    rho_xy[ind],
                    range_xy,
                    T_lsq,
                    case,
                    legend_pos_x[legend_pos],
                    legend_pos_y[legend_pos],
                )
                print(
                    (
                        "    Intercept age = %s +%s %s %s [%d%% conf.]"
                        % (
                            format(Tsi / age_unit, dignum),
                            format((Tmax - Tsi) / age_unit, dignum),
                            format((Tmin - Tsi) / age_unit, dignum),
                            age_unit_name,
                            ia_cr * 100,
                        )
                    )
                )

            if (concordia_ia_case_tw == 1) or (concordia_ia_case_tw == 2):
                legend_pos += 1
                case = 1
                Tsi, Tmax, Tmin = plot_concordia_intercept_age(
                    ax,
                    axn,
                    "tw",
                    x[ind],
                    y[ind],
                    sigma_x[ind],
                    sigma_y[ind],
                    ia_cr,
                    rho_xy[ind],
                    range_xy,
                    T_lsq,
                    case,
                    legend_pos_x[legend_pos],
                    legend_pos_y[legend_pos],
                )
                print(
                    (
                        "    Intercept age = %s +%s %s %s [%d%% conf.]"
                        % (
                            format(Tsi / age_unit, dignum),
                            format((Tmax - Tsi) / age_unit, dignum),
                            format((Tmin - Tsi) / age_unit, dignum),
                            age_unit_name,
                            ia_cr * 100,
                        )
                    )
                )

        # Sample number
        legend_data_number(
            ax, axn, legend_pos_x[legend_pos], legend_pos_y[legend_pos], ind
        )

    # ------------------------------------------------
    # C: Bar plot of 206Pb/238U ages with 1D weighted mean

    if plot_diagrams[2] == 1:
        if np.sum(plot_diagrams[0:2]) == 2:
            axn = 2
            axn_title = "c"
        elif np.sum(plot_diagrams[0:2]) == 1:
            axn = 1
            axn_title = "b"
        else:
            axn = 0
            axn_title = "a"

        print(("%s: One-dimensional bar plot") % axn_title)

        ax[axn].set_title(axn_title, loc="left", fontsize=legend_font_size + 6)
        ax[axn].set_xlim([0, N + 1])
        ax[axn].set_ylim(range_oneD_y[0], range_oneD_y[1])
        ax[axn].set_xlabel("Number of samples", fontsize=legend_font_size + 4)

        Tall, s1, label_selected = select_age_type(oneD_age_type)
        ax[axn].set_ylabel(label_selected, fontsize=legend_font_size + 4)

        legend_pos_x, legend_pos_y = calc_legend_pos(
            [[(N + 1) * 0.05, (N + 1) * 0.05], range_oneD_y]
        )
        T_owm, S_owm, MSWD_owm = plot_oneD_weighted_mean(
            ax[axn],
            oneD_age_type,
            Tall,
            s1,
            ind,
            outd,
            outd_disc,
            oneD_cr,
            legend_pos_x,
            legend_pos_y,
        )
        print(
            u"    1D weighted mean age = %s ± %s %s [%d%% conf.] (MSDW=%s)"
            % (
                format(T_owm, dignum),
                format(S_owm, dignum),
                age_unit_name,
                oneD_cr * 100,
                format(MSWD_owm, dignum),
            )
        )

        # reduced chi-squared (Spencer2016gf)
        chi2_red, res_chi2_red = calc_chi2_red(Tall[ind], s1[ind], T_owm, len(ind), opt=1)
        print(
            u"    Reduced Chi-squared = %s (%s)"
            % (format(chi2_red, dignum), res_chi2_red)
        )

    # ------------------------------------------------
    # Histogram
    # Th/U and age histogram plots

    if plot_diagrams[3] == 1:
        if np.sum(plot_diagrams[0:3]) == 3:
            axn = 3
            axn_title = "d"
        elif np.sum(plot_diagrams[0:3]) == 2:
            axn = 2
            axn_title = "c"
        elif np.sum(plot_diagrams[0:3]) == 1:
            axn = 1
            axn_title = "b"
        else:
            axn = 0
            axn_title = "a"

        print(("%s: Histogram") % axn_title)

        ax[axn].set_title(axn_title, loc="left", fontsize=legend_font_size + 6)
        ax[axn].set_xlim(range_hist_x[0], range_hist_x[1])
        Tall, s1, label_selected = select_age_type(hist_age_type)
        ax[axn].set_xlabel(label_selected, fontsize=legend_font_size + 4)

        # Optional: Th/U ratio
        if opt_Th_U:
            range_hist_y2 = [0.0, np.ceil(np.max(Th_U)) + 0.1]
            axb = ax[axn].twinx()
            plot_Th_U(
                axb,
                Th_U,
                Th_U_e,
                Tall,
                s1,
                ind,
                outd,
                outd_disc,
                oneD_cr,
                range_hist_x,
                range_hist_y2,
            )

        if opt_kde:
            plot_kde(ax[axn], range_hist_x, Tall, ind)

        plot_hist(ax[axn], Tall, ind, outd_disc, outd)

    print("All done.")

    # ------------------------------------------------
    # Output pdf file

    if "pdf" in options.driver:
        print("Saving %s" % outfile)
        plt.savefig(outfile)
    else:
        plt.show()
