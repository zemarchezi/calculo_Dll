#!/usr/bin/python2
# coding=utf-8
# author: Jose P. Marchezi
import numpy as np
#import matplotlib
import matplotlib.pyplot as plt
import os
import datetime
import matplotlib
import matplotlib.dates as mdates
from marxeDownloader import MarxeDownloader
from analysisFunc import dll_calc
#import matplotlib.ticker as ticker

# Extract single variable for .cdf file
# def cdf_singleVar(filename,variable):
#     data = pycdf.CDF(filename)
#     var = data[variable][...]
#     return var

potsdam = 1 # if 1, download the kp data from potsdam, if 0 you must import manually from spidr
# 
# timeini = datetime.datetime(2014,9,21,0,0)
# timeend = datetime.datetime(2014,9,24,12,0)



# the filename is used only when the data is from SPIDR repository
kpFilename = 'spidr_1488827144326.txt'
####
matplotlib.rc('text', usetex=True) # setting to allow LaTex comands and symbols
matplotlib.rc('font', family='serif', size=12)

# directory of the data
path = os.getcwd() #'/home/jose/MEGA/Doutorado/Progs/plot_ULF/dados/'
dataDownlDir = path + '/dados/'
plotDir = path + '/plots/'


if potsdam:
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
####
## plots
fig = plt.figure(figsize=(24, 6))
pl = fig.add_subplot(111)
sp = pl.pcolormesh(date_list2,lL,np.log10(dllt),shading='gouraud', cmap='jet') # viridis  jet
pl.set_ylabel('L star [$R_E$]')
# pl.axes.get_xaxis().set_visible(False)
pl.set_xlim(timeini, timeend)
pl.set_ylim(2, 7)
pl.set_title('Total $D_{LL}$')
pl.xaxis.set_major_formatter(mdates.DateFormatter('%Y/%m/%d'))
pl.xaxis.set_major_locator(mdates.DayLocator())
pl.xaxis.set_minor_locator(mdates.HourLocator(np.arange(0, 25, 1)))
pl.xaxis.set_tick_params(which='major', length=10, width=2, color='k')
pl.xaxis.set_tick_params(which='minor', length=5, width=2, color='k')
cbar_coord = [0.94,0.1,0.01,0.76]
cbar_ax = fig.add_axes(cbar_coord)
# cbar = fig.colorbar(sp, cax=cbar_ax, boundaries=np.linspace(0,4,64), ticks = [0,1,2,3,4], label = '$log_{10}$ (Mean Wave power) $[(nT)^2/Hz)]$')
cbar = fig.colorbar(sp, cax=cbar_ax,boundaries=np.linspace(-5,0.8,64), ticks = [-5,-4,-3,-2,-1,0], label = '$log_{10}$ ($D_{LL}$)')
cbar.set_clim(-5,0)
ax2 = pl.twinx()
ax2.step(date_list2, mkp, where='mid', color='k',linewidth=2)
# ax2.axvline(x=datetime.datetime(2014,9,22,3,0),linestyle='--', linewidth=4, color='r')
ax2.set_ylabel('Kp', color='k')
ax2.set_ylim(-5,0)
ax2.yaxis.set_ticklabels([5,4,3,2,1,0])
ax2.tick_params(axis='y', colors='red')
ax2.yaxis.label.set_color('red')
pl.set_xlim(timeini, timeend)

plt.savefig(plotDir+'DLL_Kp_%04d%02d%02d_to_%04d%02d%02d.png' % (timeini.year, timeini.month, timeini.day, timeend.year, timeend.month, timeend.day))


# plot Dll E and B in the same figure

fig = plt.figure(figsize=(24, 13))
pl = fig.add_subplot(211)
sp = pl.pcolormesh(date_list2,lL,np.log10(dllE),shading='gouraud', cmap='jet') # viridis  jet
pl.set_ylabel('L star [$R_E$]')
# pl.axes.get_xaxis().set_visible(False)
# pl.set_xlim(timeini, timeend)
pl.set_ylim(2, 7)
pl.set_title('$D_{LL}^{E}$')
pl.xaxis.set_major_formatter(mdates.DateFormatter('%d/%m'))
pl.xaxis.set_major_locator(mdates.DayLocator())
pl.xaxis.set_minor_locator(mdates.HourLocator(np.arange(0, 25, 1)))
pl.xaxis.set_tick_params(which='major', length=10, width=2, color='k')
pl.xaxis.set_tick_params(which='minor', length=5, width=2, color='k')
cbar_coord = [0.94,0.53,0.01,0.35]
cbar_ax = fig.add_axes(cbar_coord)
# cbar = fig.colorbar(sp, cax=cbar_ax, boundaries=np.linspace(0,4,64), ticks = [0,1,2,3,4], label = '$log_{10}$ (Mean Wave power) $[(nT)^2/Hz)]$')
cbar = fig.colorbar(sp, cax=cbar_ax,boundaries=np.linspace(-5,0.8,64), ticks = [-5,-4,-3,-2,-1,0], label = '$log_{10}$ ($D_{LL}$)')
cbar.set_clim(-5,0)
ax2 = pl.twinx()
ax2.step(date_list2, mkp, where='mid', color='k',linewidth=2)
ax2.axvline(x=datetime.datetime(2014,9,22,3,0),linestyle='--', linewidth=4, color='r')
ax2.set_ylabel('Kp', color='k')
ax2.set_ylim(-5,0)
ax2.yaxis.set_ticklabels([5,4,3,2,1,0])
ax2.tick_params(axis='y', colors='red')
ax2.yaxis.label.set_color('red')
pl.set_xlim(timeini, timeend)


pl = fig.add_subplot(212)
sp = pl.pcolormesh(date_list2,lL,np.log10(dllB),shading='gouraud', cmap='jet') # viridis  jet
pl.set_ylabel('L star [$R_E$]')
# pl.axes.get_xaxis().set_visible(False)
# pl.set_xlim(timeini, timeend)
pl.set_ylim(2, 7)
pl.set_title('$D_{LL}^{B}$')
pl.xaxis.set_major_formatter(mdates.DateFormatter('%d/%m'))
pl.xaxis.set_major_locator(mdates.DayLocator())
pl.xaxis.set_minor_locator(mdates.HourLocator(np.arange(0, 25, 1)))
pl.xaxis.set_tick_params(which='major', length=10, width=2, color='k')
pl.xaxis.set_tick_params(which='minor', length=5, width=2, color='k')
cbar_coord = [0.94,0.1,0.01,0.36]
cbar_ax = fig.add_axes(cbar_coord)
# cbar = fig.colorbar(sp, cax=cbar_ax, boundaries=np.linspace(0,4,64), ticks = [0,1,2,3,4], label = '$log_{10}$ (Mean Wave power) $[(nT)^2/Hz)]$')
cbar = fig.colorbar(sp, cax=cbar_ax,boundaries=np.linspace(-5,0.8,64), ticks = [-5,-4,-3,-2,-1,0], label = '$log_{10}$ ($D_{LL}$)')
cbar.set_clim(-5,0)
ax2 = pl.twinx()
ax2.step(date_list2, mkp, where='mid', color='k',linewidth=2)
ax2.axvline(x=datetime.datetime(2014,9,22,3,0),linestyle='--', linewidth=4, color='r')
ax2.set_ylabel('Kp', color='k')
ax2.set_ylim(-5,0)
ax2.yaxis.set_ticklabels([5,4,3,2,1,0])
ax2.tick_params(axis='y', colors='red')
ax2.yaxis.label.set_color('red')
pl.set_xlim(timeini, timeend)


# plt.show()

plt.savefig(plotDir+'DllB_E.png', format='png')
