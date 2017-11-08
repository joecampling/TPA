# -*- coding: utf-8 -*-
"""
Created on Sat Jul  8 11:02:12 2017

@author: Compaq
"""
from numpy import arange, exp, flipud, zeros, around
import csv
save = 1
lux = 299792458
"""Select centre frequency (um)"""
cent = 1.575
lam = cent / 1e6
cfreq = lux / lam / 1e12    # centre frequency in THz

nv = 2 ** 13    # number of time steps
dt = 4.2e-3     # time step in ps
tm = arange(nv) * dt - nv / 2 * dt
dfm = 1 / (dt * nv)
fm = arange(nv) * dfm - nv / 2 * dfm

"""set up wavelength limits (um)"""
minw = 1
maxw = 2
"""calculate respective frequencies (THz)"""
minf = lux / maxw / 1e6
maxf = lux / minw / 1e6
"""find corresponding nearest points in fmr"""
fmr = fm + cfreq
mindif1 = 1
mindif2 = 1
for a in range(nv):
    dif1 = fmr[a] - minf
    dif2 = fmr[a] - maxf
    if abs(dif1) < abs(mindif1):
        mindif1 = dif1
        minfpos = a
    if abs(dif2) < abs(mindif2):
        mindif2 = dif2
        maxfpos = a
"""work out minimum and maximum wavelength steps"""
gapswl = lux / fmr[maxfpos - 1] / 1e6 - lux / fmr[maxfpos] / 1e6
gaplwl = lux / fmr[minfpos] / 1e6 - lux / fmr[minfpos + 1] / 1e6
dlam = gapswl # wavelength step (minimum, i.e. shorter wavelength gap)
wavesteps = int((lux / fmr[minfpos] / 1e6 - lux / fmr[maxfpos] / 1e6) / dlam)
rev = flipud(fmr)
waveaxis = arange(wavesteps) * dlam + minw
eqfm = lux / waveaxis / 1e6
eqfmpos = zeros(wavesteps)
mindif = 1
for a in range(wavesteps):
    for b in range(nv):
        dif = eqfm[a] - rev[b]  # compare freq in eqfm and original freq axis
        if dif < mindif:
            rec = b
    eqfm[a] = rev[rec]  # set the eqfm to match the nearest freq step
    eqfmpos[a] = rec
    mindif = 1 

waxis = zeros((wavesteps, 2))
waxis[:,0] = waveaxis
waxis[:,1] = eqfmpos
fileopen = 'waveaxis{0}_{1}.csv'
filename = str.format(fileopen, int(around(cent * 100, 0)), wavesteps)
if save:
    file = open(filename, 'w', newline = '')
    writer = csv.writer(file)
    for row in waxis:
        writer.writerow(row)
        
    file.close()