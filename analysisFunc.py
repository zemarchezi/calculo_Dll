#!/usr/bin/python2
# coding=utf-8
# author: zemarchezi
from scipy.signal import butter, lfilter
from scipy import interpolate as interp
import numpy as np
from spacepy import pycdf
import os
from marxeDownloader import MarxeDownloader
import datetime

# Butterword filter coefficients
def butter_bandpass(lowcut, highcut, fs, order=5):
     nyq = 0.5 * fs
     low = lowcut / nyq
     high = highcut / nyq
     b, a = butter(order, [low, high], btype='band')
     return b, a

# Band-pass butterword filter
def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
     b, a = butter_bandpass(lowcut, highcut, fs, order=order)
     y = lfilter(b, a, data)
     return y

# fill the gaps with nans
def fill_nan(A):
     '''
     interpolate to fill nan values
     '''
     if np.isnan(A[0]):
         A[0] = 0.5
     inds = np.arange(A.shape[0])
     good = np.where(np.isfinite(A))
     f = interp.interp1d(inds[good], A[good], bounds_error=False)
     B = np.where(np.isfinite(A), A, f(inds))
     return B

# convert the coordinate system from gse to field aligned system
def rotate_field_fac(bx, by, bz, x, y, z, v1x, v1y, v1z):
     r = [x, y, z] / np.sqrt(x * x + y * y + z * z)
     ep = [bx, by, bz] / np.sqrt(bx * bx + by * by + bz * bz)
     ea = np.cross(ep, r) / np.linalg.norm(np.cross(ep, r))
     er = np.cross(ea, ep) / np.linalg.norm(np.cross(ea, ep))
     # apply rotation for B
     bp = (ep[0] * bx) + (ep[1] * by) + (ep[2] * bz)
     ba = (ea[0] * bx) + (ea[1] * by) + (ea[2] * bz)
     br = (er[0] * bx) + (er[1] * by) + (er[2] * bz)
     # Apply the rotation for vector V1
     v1p = (ep[0] * v1x) + (ep[1] * v1y) + (ep[2] * v1z)
     v1a = (ea[0] * v1x) + (ea[1] * v1y) + (ea[2] * v1z)
     v1r = (er[0] * v1x) + (er[1] * v1y) + (er[2] * v1z)
     # testing whether the rotation is correct
     b_fac = np.linalg.norm([bp, ba, br])
     b_orig = np.linalg.norm([bx, by, bz])
     return(bp, ba, br, v1p, v1a, v1r, b_fac, b_orig)

# read rept files

def extract_rept(files_rept):
    var_rept = ['Epoch', 'L', 'L_star', 'FEDU_Energy', 'FEDU', 'MLT']
    rep = []
    for v in var_rept:
        rep.append([])
        for l in files_rept:
            data_rept = pycdf.CDF(l)
            rep[len(rep)-1].extend(data_rept[v][...])

    for r in rep[4]:
        r[r == -9999999999999999635896294965248.000] = 'NaN'

    for i in range(0,len(rep[4])):
        temp = rep[4][i]
        temp2 = []
        for l in range(temp.shape[1]):
            temp2.append(np.nanmean(temp[:,l]))
        rep[4][i] = temp2
    for l in files_rept:
        os.remove(l)
    return rep


# Extract the variables of interest from emfsis .cdf file
def extract_emfisis(files_mag):
    var_mag = ['Epoch', 'Mag', 'coordinates']
    emf = []
    for v in var_mag:
        emf.append([])
        for l in files_mag:
            data_mag = pycdf.CDF(l)
            emf[len(emf)-1].extend(data_mag[v][...])

    for i in range(1,len(emf)):
        temp = np.asmatrix(emf[i])
        temp[temp == -99999.898] = 'NaN'
        emf[i] = temp

    return emf

# Extract the variables of interest from efw .cdf file
def extract_efw(files_ele):
    var_ele = ['epoch', 'efield_inertial_frame_mgse', 'spinaxis_gse']
    efw = []
    for v in var_ele:
        efw.append([])
        for l in files_ele:
            data_ele = pycdf.CDF(l)
            efw[len(efw)-1].extend(data_ele[v][...])

    for i in range(1,len(efw)):
        temp = np.asmatrix(efw[i])
        temp[temp == -9999999848243207295109594873856.000] = 'NaN'
        efw[i] = temp
    return efw
##
##
##
## DLL Calculation
def dll_calc(timeini, timeend, downloadData, filename = 'none'):
    potsdam = downloadData # if 1, download the kp data from potsdam, if 0 you must import manually from spidr
    timeini = timeini
    timeend = timeend

    # the filename is used only when the data is from SPIDR repository
    ###

    # directory of the data
    path = os.getcwd() #'/home/jose/MEGA/Doutorado/Progs/plot_ULF/dados/'
    dataDownlDir = path + '/dados/'
    plotDir = path + '/plots/'


    if downloadData == 1:
        # Download data
        str_temp = 'kp%s%02d.tab' % (str(timeini.year)[2:], timeini.month)
        # define the directory and host in the ftp
        host = 'ftp.gfz-potsdam.de'
        working_directory = '/pub/home/obs/kp-ap/tab/'
        #### download of the data #################
        mx = MarxeDownloader(host)
        # Connecting
        mx.connect()
        # Set downloaded data directory
        mx.set_output_directory(dataDownlDir)
        # Set FTP directory to download
        mx.set_directory(working_directory)
        # Download single data
        mx.download_one_data(str_temp)

        f = open(dataDownlDir + str_temp, 'r')
        data = f.readlines()
        f.close()
        data = data[:-4]
        useKp = data[(timeini.day-1):(timeend.day)]
        usekp2 = []
        for i in range(0, len(useKp)):
            temp = useKp[i].split('  ')
            usekp2.append([])
            usekp2[i].extend(temp[1].split(' '))
            usekp2[i].extend(temp[2].split(' '))
            totKpDays = []
        for ii in usekp2:
            totKpDays.extend(ii)
        for i in range(0, len(totKpDays)):
            temp = totKpDays[i]
            if temp[1] == 'o':
                totKpDays[i] = float(totKpDays[i][0])
            if temp[1] == '-':
                totKpDays[i] = float(totKpDays[i][0]) - 1./3
            if temp[1] == '+':
                totKpDays[i] = float(totKpDays[i][0]) + 1./3
    else:
        kpFilename = filename
        kp = np.loadtxt(dataDownlDir+kpFilename,delimiter=' ', skiprows=16, usecols=(2,2), unpack=True)
        totKpDays = kp[0]

    kp = totKpDays
    mkp = []
    # invert the signal of kp in order to plot togheter with the dll
    for i in totKpDays:
        mkp.append(i*-1)

    ttt = np.linspace(1,4,len(kp))

    # create a date array
    ran = (len(kp) / 8) * 24
    date_list2 = [datetime.datetime(timeini.year, timeini.month, timeini.day, 0, 0) + datetime.timedelta(hours=x) for x in range(0,ran,3)]


    # perfor the Dll calculation
    lL = np.linspace(2,7,100)
    dllE = []
    dllB = []
    dllt = []
    for l in lL:
        tempDllB = []
        tempDllE = []
        tempDllt = []
        for k in kp:
            expB = (-0.0327 * np.power(l,2)) + 0.625*l - 0.0108 * np.power(k,2) + 0.499 * k
            expE = (0.217 * l ) + (0.461 * k)
            tempB = 6.62e-13 * np.power(l,8) * np.power(10, expB)
            tempE = 2.16e-8 * l**6 * 10**(expE)
            tempt = tempB + tempE
            tempDllB.append(tempB)
            tempDllE.append(tempE)
            tempDllt.append(tempt)
        dllB.append(tempDllB)
        dllE.append(tempDllE)
        dllt.append(tempDllt)

    return(dllB,dllE,dllt, date_list2, lL, kp)
