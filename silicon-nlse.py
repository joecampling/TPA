# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 15:09:22 2017

@author: Joseph Campling, Optoelectronics Reseach Centre (ORC), University of Southampton
"""

from numpy import sqrt, pi, zeros, ones, append, arange, exp, log10, sin, log, abs, seterr, real, max, copy
from numpy.fft import ifftshift as fftshift, fftshift as ifftshift
from numpy.random import randn, poisson
from scipy.fftpack import ifft as fft, fft as ifft
from matplotlib.pyplot import plot, show, xlabel, ylabel, figure, close, xlim, ylim, tick_params, imshow \
    , subplot2grid, colorbar, legend
from matplotlib.ticker import FormatStrFormatter
from scipy.integrate import simps
from math import factorial
from nlse_functions import progress, fca, tpadis, sech, SSFMTR2, third, fca4
import csv
import glob
seterr(invalid = 'ignore')
seterr(divide = 'ignore')
#sys.exit()
compare = 1
include3pa = 1
linear_loss = 0
study2only = 0
posbeta2 = 0
n2mult = 1  # factor to increase n2 by
linecol = 'k'
noiselevel = 0.5e-5
d3pa_2 = 0
d3pa_3 = 0
Cj3 = 0 
fRsi = 0.043#0.03409082702482847  # fractional contribution of Raman in silicon

if include3pa:
    t3pa = 0.0018e-24
    d3pa = -3e-30##-7.5788e-29 
    d3pa_2 = 0#-2e-33 #6.65e-33
    d3pa_3 = 0#6e-35 #6e-36
else:
    t3pa = 0
    d3pa = 0
    d3pa_2 = 0
    d3pa_3 = 0

"""Choose Sellmeier equation"""
Sel = 1

"""Close all figures?"""
closeall = 0

"""Which figures to show?"""
figstoshow = {'Figure1': 1, 'Figure2': 1, 'Figure3': 0, 'Figure4': 0, 'Figure5': 0}

M = 2 ** 8               # number of steps in z

"""Select material (1 = a-Si-H, 2 = p-Si, default = a-Si-H)"""
material = 1

"""Input carrier frequency in um"""
carr = 1.575

"""Pulse width"""
tau = 0.3 / 1.665   # T0 (ps)

"""Input power of modes"""
Average_power = 20 # mW
p0 = 75
p1 = 0
p2 = 0
p3 = 0

"""Length of fiber (mm)"""
length = 4

"""Fiber diameter (nm)"""
D = 1700   

att_dB_over_cm = 0.8 #+ (2.4 - carr) / 0.85 * 2

"""Test soliton?"""
soliton = 0

matname = 'a-Si-H' # default
if material == 2:
    matname = 'p-Si'
"""
Import n2 values for a-Si-H
n2dic = {}
file = open('n2_Si.csv')
reader = csv.reader(file)
for row in reader:
    n2dic[eval(row[0])] = eval(row[1])
file.close()

if carr == 2.15:
    n2 = 1.388e-17 
elif carr == 3.7:
    n2 = 6e-18
else:
    n2 = n2dic[carr]
if material == 2: # reduce n2 values for p-Si
    n2 /= 2.1 # / 2.1 to match avg. chi3 curve in Jalali at al. / 1.383 to get Bristow curve
"""
if closeall:
    close('all')

tpad = 0
TD = 1
if not(figstoshow['Figure3'] or figstoshow['Figure4'] or figstoshow['Figure5']):
    TD = 0  # save as a 2D grid over all of z

"""Load in effective areas
openfile = 'Areas_{0}_{1}.csv'
filename = str.format(openfile, matname, D)
file = open(filename)
reader = csv.reader(file)
Areas = {}
for row in reader:
    Areas[eval(row[0])] = eval(row[1])

# load in betas
beta2s = {}
openfile = 'beta2s_{0}_{1}.csv'
filename = str.format(openfile, matname, D)
file = open(filename)
reader = csv.reader(file)
for row in reader:
    beta2s[eval(row[0])] = eval(row[1])
file.close()

beta3s = {}
openfile = 'beta3s_{0}_{1}.csv'
filename = str.format(openfile, matname, D)
file = open(filename)
reader = csv.reader(file)
for row in reader:
    beta3s[eval(row[0])] = eval(row[1])
file.close()
"""
I12 = 144   # power for intensity = 1.2GW/cm**2
I60 = 720   # 6.0GW/cm**2
I125 = 1500 # 12.5GW/cm**2

tpa = 0.2e-12
t4pa = 0#3.6e-41
fa = 1     # turns fca on/off
fd = 1     # turns fcd on/off
dtpa = 2.378e-14 / 2
# fixed parameters
lux = 299792458             # speed of light
h = 6.62607004e-34          # Planck's constant
fR = 0.18                   # fractional Raman contribution in silica
fRn = 0.245                 # fraction of revised Raman response


"""Grid setup"""
lz = length / 1e3
dt = 4.2e-3 # 2.8e-3                 # time step [ps]
dts = dt / 1e12                # time step [s] - for integrating over pulse
nv = 2 ** 13                # number of pixels in time
if M < 2 ** 12:
    savepoints = M
else:
    savepoints = 2 ** 12              # only if nv * M <= 2^26
dz = lz / M                 # space step [m]**-1
z = (arange(M) + 1) * dz    # space axis
dfm = 1 / (nv * dt)         # frequency step [Thz]
fm = (arange(nv) - nv / 2) * dfm  # frequency axis
wm = 2 * pi * fm            # angular frequency axis (rad THz)
wmm = wm * 1e12             # angular frequency axis (rad Hz)
lam = carr * 1e-6           # centre wavelength
cfreq = lux / lam / 1e12    # centre frequency (THz)
wave = lux / (cfreq + fm) / 1e3     # wavelength axis (nm)

# Turn Raman and self-steepening on/off
Ram = 1     # turn Raman on / off
SSt = 1    # turn self-steepening on/off

"""Mode overlaps"""
R11=1 / 1.4e-12 # not used, naive value assuming zero z-component of modal field
R22=1 / 1.88257381982e-12
R33=1 / 9.58984790412e-12
R44=1 / 9.66104512147e-12
R12=1 / 2.31936196395e-12
R13=1 / 15.016931334399459e-12
R14=1 / 15.106457720226712e-12
R23=1 / 12.571752185265419e-12
R24=1 / 13.04120857793736e-12
R34=1 / 12.274110915706891e-12

t0 = tau * 1e-12    # T0 (s)
k0 = 2 * pi / lam
sigma = 1e-20 * carr / 1.55
kc = 1.35e-27
muf = 2 * kc * k0 / sigma
om0 = 2 * pi * cfreq
nu0 = lux / lam

# betas
# D = 4e-12 / (1e-9 * 1e3)  # dispersion in ps/nm-km
# beta2 = -lam ** 2 / (2 * pi * lux) * D * 1e24  # Beta2 in ps**2/m
beta_dif_12 = (1.4142002362758894e-08 - 1.2518730802477717e-08) * 1e12
beta_dif_13 = 2.53505e-10 *1e12
beta_dif_14 = 5.07796e-10 *1e12

beta4 = 0#-2.6671732822827566e-5 #-4.081352915829911e-5
beta5 = 0#7.740634783138472e-7
beta6 = 0#-1.130376215162516e-8
beta7 = 0#1.262298363451786e-10
beta8 = 0#-3.081409500391959e-13
beta9 = 0#-4.489874605017346e-14
beta10 = 0#2.000814953490767e-15
beta11 = 0#-5.440777102870384e-17
beta12 = 0#4.348011401440699e-19
beta13 = 0#8.759189864679899e-20
beta14 = 0#-6.108426498790149e-21
beta15 = 0#-5.899364193804093e-23
beta16 = 0#-5.074522239137319e-23
beta17 = 0#1.196528342487147e-24

if carr < 3:
    R11 = 1 / Areas[carr]
    beta2 = beta2s[carr] * 1e24
if carr == 2.15:
    beta2 = -5.864960690387303e-1
    beta3 = 1.935130836756176e-3 
    beta4 = -8.672265297710004e-6
    beta5 = 4.939086331219052e-8
    beta6 = -3.423605644999189e-10
    beta7 = 2.800567451951438e-12
elif carr == 3.7:
    R11 = 1 / 1.15e-12
    beta2 = -8.768377382573194e-2
    beta3 = 2.040986425170804e-3
else:
    beta3 = beta3s[carr] * 1e36

Rs = [[R11, R12, R13, R14], [R12, R22, R23, R24], [R13, R23, R33, R34], [R14, R24, R34, R44]] # matrix of mode overlaps
dAf = 2e-5 / R11     # dAeff / dw
beta22 = 0
beta23 = 0
beta24 = 0.1

absor = 1 / 2 * (att_dB_over_cm / 10 * log(10) / 1e-2)

"""Set up array for wavelength-dependent linear loss"""
alpha = zeros(nv)
fm_2 = lux / 2e-6 / 1e12
fm_24 = lux / 2.4e-6 / 1e12
fm_26 = lux / 2.6e-6 / 1e12
no_steps1 = (fm_24 - fm_26) / dfm
no_steps = (fm_2 - fm_24) / dfm
dB_step0 = 0.5 / no_steps1
#if section:
 #   dB_step0 *= 40
no_steps_mid = (fm_2 - fm_24) / dfm
dB_step_mid = 3.5 / no_steps_mid #3.5 / no_steps
dB_step2 = 3.6 / 2476
no_steps_max = (fm_24 - cfreq - fm[0]) / dfm
max_dB = no_steps_max * dB_step0

for a in range(nv):
    if lux / (cfreq + fm[a]) / 1e3 > 2800 and lux / (cfreq + fm[a]) / 1e3 < 2800:
        alpha[a] = 100
    else:
        alpha[a] = 0.8
        """
        alpha[a] = max_dB - a * dB_step0 + 4.5
        if a > no_steps_max:
            alpha[a] = 4.5 - (a - no_steps_max) * dB_step_mid
        if a > (no_steps_max + no_steps_mid):
            alpha[a] = 1 + (a - no_steps_max - no_steps_mid) * dB_step2
        """
if linear_loss:
    absor = 1 / 2 * (alpha / 10 * log(10) / 1e-2)
    
tm = (arange(nv) - nv / 2) * dt  # symmetrical time axis
tm0 = tm + nv / 2 * dt # time axis starting at 0
tmh = append(zeros(int(nv / 2)), ones(int(nv / 2)))  # heaviside step function as array
n2 *= n2mult
dn2 =  2.564e-32 / n2 # dn2 / dom
if carr < 2.3:
    dn2 = 3.083e-32 / n2 # higher dn2/dom for wavelengths < 2.3um
ar = 1 / R11 # effective area of fundamental mode
gam = (k0 * n2 + 1j /2 * tpa) # gamma' (includes tpa but not Aeff)
regam = real(gam) * R11 # real part of gamma (incl, Aeff, i.e. standard non-semiconductor-gamma)
DW = lux / ((-3 * beta2 / beta3 + regam * p0 * beta3 / 3 / beta2 ** 2) / 2 / pi + cfreq) / 1e6 # wavelength of DW
ld = tau ** 2 / abs(beta2) # dispersion length
N = sqrt(regam * p0 * t0 ** 2 / abs(beta2 * 1e-24))     # soliton order
sfl = ld / N    # soliton fission length
amp = sqrt(p0) # amplitude of peak input fundamental mode
lnl = 1 / real(gam) / amp  # nonlinear length
I0 = p0 * R11 # intensity of each mode
I1 = p1 * R22
I2 = p2 * R33
I3 = p3 * R44
cond = 3 * h * nu0 / sigma / t0     # condition to compare with I0. If I0 << cond, FCA is negligible
sqr = append(zeros(int(7 * nv / 16)), ones(int(nv / 8)))
sqr = append(sqr, zeros(int(7 * nv / 16)))
sqrf = ifft(sqr)
Cc = tpa / 2 / h / nu0
Cj = Cc * sigma * (fa + fd * 1j * muf)
sigma4 = 8.26e-21
mu4 = 3.16
C4 = t4pa / 4 / h / nu0
Cj4 = C4 * sigma4 * (fa + fd * 1j * mu4)

if TD == 1:
    vv = zeros((savepoints, nv), dtype=complex)  # matrix to save the data in time and space
    vvF = zeros((savepoints, nv), dtype=complex)  # matrix to save the data in frequency and space

C = 0.78135485849231634 #1.2 chirp parameter
vs = zeros((4,nv), dtype = complex) # initialise arrays for pulses
vf = zeros((4,nv), dtype = complex)
vp1 = zeros((4,nv), dtype = complex)
noise = randn(nv) * 1e-5 # ensures noise floor remains constant throughout propagation

"""Set up pulses"""
vs[0] = sqrt(p0) * exp(- tm ** 2 / (2 * tau ** 2)) * exp(-1j * C * tm ** 2 / (tau ** 2))
vs[1] = sqrt(p1) * exp(- tm ** 2 / (2 * tau ** 2)) * exp(-1j * C * tm ** 2 / tau ** 2)
vs[2] = sqrt(p2) * sech(tm / tau) * exp(-1j * C * tm ** 2 / tau ** 2)
vs[3] = sqrt(p3) * sech(tm / tau) * exp(-1j * C * tm ** 2 / tau ** 2)

"""soliton parameters for testing"""
if soliton:
    beta3 = 0
    #absor = 0
    beta2 = 1
    p0 = abs(beta2 / 1e24) / regam / (t0) ** 2
    gam = real(gam)
    dn2 = 0
    dAf = 0
    amp = sqrt(p0)
    vs = zeros((4,nv), dtype = complex)
    vs[0] = amp * sech(tm / tau)
    sigmaf = zeros((4, nv), dtype = complex)
    fRsi = 0
    Cj = 0 + 0j
    
vp1[:, :] = vs[:, :] # snapshot of initial pulses in time
noisefreq = poisson(1, nv) * noiselevel + 1j * poisson(1, nv) * noiselevel
vf1 = fftshift(fft(ifftshift(vp1)))
vf1[0] += noisefreq # inital pulses in frequency
vp1 = ifftshift(ifft(fftshift(vf1)))
vs[:, :] = vp1[:, :]
init = simps(abs(vp1[0]) ** 2, None, dt)
Nct = fca(Cc, vs, Rs, dts) # test of fca function
Nct4 = fca4(Cj4, vs, Rs, dts)

"""Raman response function parameters"""
tau1 = 12.2e-3  # values for Raman shift in silica (fs)
tau2 = 32e-3
taub = 96e-3  # values for upgraded Raman response function
fb = 0.21
hR = (tau1 ** (-2) + tau2 ** (-2)) * tau1 * exp(-(tm * tmh) / tau2) * sin((tm * tmh) / tau1)
hRn = ((1 - fb) * hR + fb * ((2 * taub - (tm * tmh)) / taub ** 2) * exp(-(tm * tmh) / taub)) * tmh
tR1 = 10e-3
tR2 = 3
hRsi = 1 / tR1 * exp(-(tm * tmh) / tR2) * sin((tm * tmh) / tR1)
hRsi0 = 1 / tR1 * exp(-(tm0) / tR2) * sin((tm0) / tR1)
TRSi = Ram * fRsi * simps(tm0 * hRsi0, None, dt) # First moment of silicon Raman response
hRsi = hRsi / simps(hRsi, None, dt)
hRn = hRn / simps(hRn, None, dt)  # normalisation of Raman response so integral = 1
hR = hR / simps(hR, None, dt)
g = ifft(hR) # ifft of Raman function
gn = ifft(hRn)
gsi = ifft(hRsi)
"""
if Sel == 2.47:
    beta2 = -5.839056652212266e-1 # -5.851047080588314e-1
    beta3 = 5.574403740343147e-3 #5.571240175556175e-3
    beta4 = 0#-2.582212512092238e-5 #-2.583013493126488e-5
    beta5 = 0#1.853917038482567e-7 #1.851658230191479e-7
    beta6 = 0#-1.619757440465080e-9 #-1.616480029750378e-9
    beta7 = 0#1.735381077049374e-11 #1.727071768106941e-11
    beta8 = 0#-2.159198455528519e-13 #-2.143576699441499e-13
    beta9 = 0#3.135886050779296e-15 #3.099058944477383e-15
    beta10 = 0#-5.149056776728474e-17 #-5.037269164130355e-17
    beta11 = 0#9.557799241922929e-19
    beta12 = 0#-1.965429245993785e-20
    beta13 = 0#4.482039320885132e-22
"""
"""Set up linear operator for each mode"""
if posbeta2:
    beta2 = abs(beta2)
if study2only:
    beta3 = 0
L = zeros((4, nv), dtype = complex)
L[0] = fftshift(1j * (beta2 / 2 * wm ** 2 + beta3 / 6 * wm ** 3 + beta4 / 24 * wm ** 4 + beta5 / 120 * wm ** 5 
     + beta6 / 720 * wm ** 6 + beta7 / 5040 * wm ** 7 + beta8 / factorial(8) * wm ** 8 + beta9 / factorial(9) * wm ** 9 \
        + beta10 / factorial(10) * wm ** 10 + beta11 / factorial(11) * wm ** 11 + beta12 / factorial(12) * wm ** 12 \
        + beta13 / factorial(13) * wm ** 13 + beta14 / factorial(14) * wm ** 14 + beta15 / factorial(15) * wm ** 15 \
        + beta16 / factorial(16) * wm ** 16 + beta17 / factorial(17) * wm ** 17) - absor)  # linear operator (shifted to ifft space)
L[1] = fftshift(1j * (beta_dif_12 * wm + (beta22 / 2 * wm ** 2)) - absor)
L[2] = fftshift(1j * (beta_dif_13 * wm + (beta23 / 2 * wm ** 2)) - absor)
L[3] = fftshift(1j * (beta_dif_14 * wm + (beta24 / 2 * wm ** 2)) - absor)

tomod = M / savepoints
fine = True
for u in range(4):
    vf[u] = fft(ifftshift(vs[u])) * exp(dz / 2 * L[u])  # First half linear step
    vs[u] = ifftshift(ifft(vf[u]))

if tpad:
    
    for num in range(M):
        vs = tpadis(vs, Rs, gsi, fRsi, gam, wm, dt, om0, dz, L, dn2, dAf, dtpa, tpa, Cj, dts)
        progress(M,num)
        
elif compare:
    for num in range(M):
        vs = third(vs, Rs, gsi, fRsi, gam, wm, dt, om0, dz, L, dn2, dAf, dtpa, tpa, Cj, dts, t3pa, d3pa, Cj3, nv)
        progress(M, num)
        if TD:
            if num % tomod == 0:
                vv[int(num / tomod)] = vs[0] # if there are more than 2^12 spatial steps, grid has to miss out intermediate steps
           
else:
    for num in range(M):
        vs, fine = SSFMTR2(vs, Rs, TRSi, gam, wm, dt, om0, dz, L, dn2, dAf, dtpa, tpa, Cj, dts, t3pa, d3pa, d3pa_2 \
                           , d3pa_3, Cj3, nv, num, p0, soliton, M, fine)
        progress(M, num)
        if TD:
            if num % tomod == 0:
                vv[int(num / tomod)] = vs[0] # if there are more than 2^12 spatial steps, grid has to miss out intermediate steps

#vs += 0.05 * vp1
vvF = ifftshift(ifft(fftshift(vv)))                            
vf = fftshift(fft(ifftshift(vs)))
vsum = vs[0] + vs[1] + vs[2] + vs[3]
vfsum = fftshift(fft(ifftshift(vsum)))
vfsum1 = vfsum
if study2only:
    vf2only = copy(vfsum)

#carr = 2.42 # cheat to use wavelength axis centred at 2.42um rather than 2.45)

if figstoshow['Figure4'] or figstoshow['Figure5']:
    """Read in frequency to wavelength axis conversion csv file (output from frequency_to_wavelength_axis.py)"""
    openfile = 'waveaxis{}_'
    filename = str.format(openfile, int(carr * 100)) # filename contains centre wavelength in um without . eg 242 for 2.42um
    globby = glob.glob('./' + filename + '*') # glob finds the csv file, which has the number of wavelength steps as the
    file = open(globby[0])                      # end of its filename eg 'waveaxis242_18948.csv'
    reader = csv.reader(file)
    jtn = globby[0].replace('.\\' + filename, '') # 2 steps to remove the rest of the filename (incl. path from glob)
    jtn = jtn.replace('.csv', '')
    wavesteps = int(jtn) # extract the number of wavelength steps from the filename
    wlarr = zeros((savepoints, wavesteps), dtype = complex)
    eqfmpos = zeros(wavesteps)
    rownum = 0
    for row in reader:
        eqfmpos[rownum] = row[1]
        rownum += 1
    file.close()
    for a in range(savepoints):
        for b in range(wavesteps):
            wlarr[a, b] = vvF[a, int(eqfmpos[b])]

"""Smooth out noise fluctuations in relative error calc"""
if soliton:
    rel = zeros(nv)
    cond = 1
    if p0 < 4:
        cond = p0 / 4
    for num in range(nv):
        if abs(vs[0][num]) >= cond and abs(vp1[0][num]) >= cond:
            rel[num] = abs(1 - abs(vs[0][num]) ** 2 / abs(vp1[0][num]) ** 2)
    sol = False
    if max(rel) < 8e-3:
        sol = True
    print("\nSoliton is ", sol, ". h = ", dz)
    print('Overall relative error is', '%.3E' % max(rel))

"""Energy conservation check"""

fin = simps(abs(vs[0]) ** 2, None, dt)
loss = (init - fin) / init * 100
print ("Energy loss = ", '%.6f' % loss, "%")

"""Print expected DW wavelength and soliton fission length"""
print('DW should appear at', DW, 'um')
print('Soliton fission length:', sfl)

"""Measure SC width @ -30dB"""
thirtydB = []
compare = zeros(3)
for a in range(nv - 2):
    for b in range(3):
        compare[b] = 10 * log10(abs(vfsum[a + b]) ** 2 / max(abs(vfsum) ** 2))
    if compare[0] > -30 and compare[2] < -30 or compare[0] < -30 and compare[2] > -30:
        thirtydB.append(lux / (fm[a + 1] + cfreq) * 1e-6)
list.reverse(thirtydB)
#print('-30dB crossing points:')
#for a in thirtydB:
 #   print('%.2f' % a)

"""Figures"""

if figstoshow['Figure1']:
    """Fig. 1: Final spectrum in wavelength"""
    """
    fig1 = figure()
    axis = fig1.add_axes([0.15, 0.15, 0.8, 0.8])
    axis.ticklabel_format(useOffset=False)
    xlim(1500, 3500)
    #ylim(-39, 4)
    xlabel('Wavelength (nm)', size = 21)
    ylabel(r'Spectral density (dB / nm)', size = 21)
    tick_params(labelsize = 19)
    axis.xaxis.set_major_formatter(FormatStrFormatter('%g'))
    axis.yaxis.set_major_formatter(FormatStrFormatter('%g'))
    plot(wave, (abs(vfsum) ** 2/ max(abs(vfsum) ** 2)), 'r', linewidth = 3)
    plot(wave, (abs(vf1[0]) ** 2 / max(abs(vfsum) ** 2)), 'g', linewidth = 3, linestyle = (0, (2, 4)))
    """
    
    fig1 = figure()
    axis = fig1.add_axes([0.15, 0.15, 0.8, 0.8])
    axis.ticklabel_format(useOffset=False)
    xlim(1500, 3500)
    ylim(-39, 4)
    xlabel('Wavelength (nm)', size = 21)
    ylabel(r'Spectral density (dB / nm)', size = 21)
    tick_params(labelsize = 19)
    axis.xaxis.set_major_formatter(FormatStrFormatter('%g'))
    axis.yaxis.set_major_formatter(FormatStrFormatter('%g'))
    
    plot(wave, 10 * log10(abs(vfsum) ** 2 / max(abs(vfsum) ** 2)), linecol, linewidth = 3)
    plot(wave, 10 * log10(abs(vf1[0]) ** 2 / max(abs(vfsum) ** 2)), 'g', linewidth = 3, linestyle = (0, (2, 4)))
    

if figstoshow['Figure2']:
    """Fig. 2: Final spectrum in frequency"""
    fig2 = figure()
    plot(fm, 10 * log10(abs(vfsum) ** 2), fm, 10 * log10(abs(vf1[0]) ** 2))
    xlabel('Frequency (THz)')
    ylabel('Spectral density (dB)')

if figstoshow['Figure3']:
    """Fig. 3: Evolution in time"""
    fig3 = figure()
    axis3 = fig3.add_axes([0.13, 0.13, 0.8, 0.82])
    test = imshow(10*log10(abs(vv)**2), extent = [min(tm), max(tm), 0, lz * 1e3], aspect = 'auto', origin = 'lower',
                  cmap = 'nipy_spectral')
    xlim(-2, 2)
    xlabel('Time (ps)', size = 21)
    ylabel('Distance (mm)', size = 21)
    cbar3 = fig3.colorbar(test, shrink = 0.5, aspect = 5)
    cbar3.set_label('Power (dBm)', size = 20, rotation = 270, labelpad = 22)
    tick_params(labelsize = 20)

if figstoshow['Figure4']:
    """Fig. 4: Evolution of spectrum in wavelength"""
    fig4 = figure()
    axis4 = fig4.add_axes([0.13, 0.13, 0.8, 0.82])
    test4 = imshow(10 * log10(abs(wlarr) ** 2 / max(abs(wlarr) ** 2)), extent = [1.4, 4.1, 0, lz * 1e3], 
                   aspect = 'auto', origin = 'lower', cmap = 'nipy_spectral')
    xlabel(r'Wavelength ($\mu$m)', size = 21)
    ylabel('Distance (mm)', size = 21)
    cbar4 = fig4.colorbar(test4, shrink = 0.5, aspect = 5)
    cbar4.set_label(r'Spectral density (dB / $\mu$m)', size = 20, rotation = 270, labelpad = 22)
    tick_params(labelsize = 20)

if figstoshow['Figure5']:
    """Fig. 5: Multiple figure output for publication"""
    fig5 = figure()
    ax1 = subplot2grid((100, 100), (0, 0), rowspan = 35, colspan = 95)
    plot(wave, 10 * log10(abs(vfsum) ** 2 / max(abs(vf1[0]) ** 2)), linecol, linewidth = 3)
    plot(wave, 10 * log10(abs(vf1[0]) ** 2 / max(abs(vf1[0]) ** 2)), 'g', linewidth = 3, linestyle = (0, (2, 4)))
    #plot(wave, abs(vfsum) ** 2 * 30 / max(abs(vfsum) ** 2) -70, 'r', linewidth = 3)
    #plot(wave, 10 * log10(abs(vf_pos) ** 2 / max(abs(vfsum) ** 2)) -36, 'b', linewidth = 3)
    ax1.set_yticks([0, -20, -40])
    ax1.set_xlim(1500, 3500)
    ax1.set_ylim(-70, 2)
    ax1.set_xlabel('Wavelength (nm)', size = 18)
    ax1.set_ylabel('Power (dB)', size = 18)
    tick_params(labelsize = 18)
    ax4 = ax1.twinx()
    ax4.plot(wave, abs(vf75W) ** 2 / max(abs(vf75W) ** 2), 'r', linewidth = 3)
    ax4.set_xlim(1500, 3600)
    ax4.set_ylim(0, 1.3)
    ax4.set_yticks([0, 0.5, 1])
    ax4.set_yticklabels(['0', '0.5', '1'], color = 'r')
    ax4.set_ylabel('Power (a.u.)', rotation = 270, size = 18, labelpad = 25, color = 'r')
    tick_params(labelsize = 18)
    ax2 = subplot2grid((100, 100), (55, 0), rowspan = 45, colspan = 40)
    tspec = imshow(10 * log10(abs(vv) ** 2 / max(abs(vv) ** 2)), extent = [min(tm), max(tm), 0, lz * 1e3], 
                  aspect = 'auto', origin = 'lower', cmap = 'nipy_spectral')
    cbart = colorbar(tspec, ax = ax2)
    cbart.set_label('Power (dB)', size = 18, rotation = 270, labelpad = 20)
    cbart.set_ticks([0, -50, -100, -150])
    xlim(-2, 2)
    xlabel('Time (ps)', size = 18)
    ylabel('Distance (mm)', size = 18)
    tick_params(labelsize = 18)
    ax3 = subplot2grid((100, 100), (55, 50), rowspan = 45, colspan = 50)
    wlspec = imshow(10 * log10(abs(wlarr) ** 2 / max(abs(wlarr) ** 2)), extent = [1.4, 4.1, 0, lz * 1e3], 
                   aspect = 'auto', origin = 'lower', cmap = 'nipy_spectral')
    cbarwl = colorbar(wlspec, ax = ax3)
    cbarwl.set_label('Power (dB)', size = 18, rotation = 270, labelpad = 20)
    cbarwl.set_ticks([0, -50, -100, -150, -200])
    xlabel(r'Wavelength ($\mu$m)', size = 18)
    tick_params(labelsize = 18, labelleft = 'off')
    ax1.set_title('(a)', loc = 'left', size = 18)
    ax2.set_title('(b)', loc = 'left', size = 18)
    ax3.set_title('(c)', loc = 'left', size = 18)
    ax2.set_xticks([-2, -1, 0, 1, 2])
    ax2.set_yticks([0, 1, 2, 3, 4])
    ax3.set_xticks([1.5, 2.5, 3.5])
    ax1.annotate('DW', xy = (1634, -47), fontsize = 17)
    ax1.annotate('', xy = (1620, -21), xytext = (1673, -39), arrowprops = dict(facecolor = 'black', width = 1 \
                 , headwidth = 10), fontsize = 17)
    ax1.annotate('',xy = (1707, -65), xytext = (1696, -49), arrowprops = dict(facecolor = 'black', width = 1 \
                 , headwidth = 10), fontsize = 17)
    ax1.annotate('SB2', xy =(2000, -29), xytext = (1954, -50), arrowprops = dict(facecolor = 'black', width = 1 \
                 , headwidth = 10), fontsize = 17)
    ax1.annotate('', xy =(2083, -65), xytext = (2044, -52), arrowprops = dict(facecolor = 'black', width = 1 \
                 , headwidth = 10), fontsize = 17)
    ax1.annotate('SB1', xy = (2730, -12), xytext = (2797, -43), arrowprops = dict(facecolor = 'black', width = 1 \
                 , headwidth = 10), fontsize = 17)
    ax1.annotate('', xy = (2765, -57), xytext = (2797, -45), arrowprops = dict(facecolor = 'black', width = 1 \
                 , headwidth = 10), fontsize = 17)
    ax1.annotate('SB1', xy = (2248, -16), xytext = (2180, -38), arrowprops = dict(facecolor = 'black', width = 1 \
                 , headwidth = 10), fontsize = 17)
    ax1.annotate('', xy = (2324, -53), xytext = (2262, -40), arrowprops = dict(facecolor = 'black', width = 1 \
                 , headwidth = 10), fontsize = 17)

"""Fig. 6: For soliton test, relative error between output and input in time domain"""
if soliton:
    fig6 = figure()
    plot(tm, rel)
    xlabel('Time (ps)')
    ylabel('Relative error')
    xlim(-1, 1)

show()
"""
figure()
plot(wave, 10 * log10(abs(vf1mW) ** 2 / max(abs(vf1mW) ** 2)), 'k', label = '1mW')
plot(wave, 10 * log10(abs(vf5mW) ** 2 / max(abs(vf5mW) ** 2)), 'b', label = '5mW')
plot(wave, 10 * log10(abs(vf10mW) ** 2 / max(abs(vf10mW) ** 2)), 'g', label = '10mW')
plot(wave, 10 * log10(abs(vf20mW) ** 2 / max(abs(vf20mW) ** 2)), 'c', label = '20mW')
plot(wave, 10 * log10(abs(vf40mW) ** 2 / max(abs(vf40mW) ** 2)), 'p', label = '40mW')
plot(wave, 10 * log10(abs(vf60mW) ** 2 / max(abs(vf60mW) ** 2)), 'r', label = '60mW')
xlim(1500, 2500)
ylim(0, -25)
legend()
"""