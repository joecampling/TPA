# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 14:38:04 2017

@author: Joseph Campling, Optoelectronics Reseach Centre (ORC), University of Southampton
IAP method by Hult (and nonlinear operators are 12 and 14 from Deiterding et al)
"""

from numpy import exp, abs, conj, pi, real, imag, shape, zeros, copy, cosh, vectorize, max
from numpy.fft import ifftshift as fftshift, fftshift as ifftshift
from scipy.fftpack import ifft as fft, fft as ifft
from scipy.integrate import cumtrapz
from numpy.random import randn
import sys
lux = 299792458             # speed of light
lam = 1550e-9               # centre wavelength
hev = 4.135667662e-15       # Plack's constant (eV)
hevb = hev / 2 / pi         # h-bar (eV)
Eig = 1.12                  # Silicon indirect bandgap (eV)
Eigm = Eig * 1.6e-19        # bandgap in mks units
K = 7.895e-66               # Curve-fitting parameter
n = 3.47752931058           # refractive index of silicon @ 1.55um
con = K / n ** 2 / Eigm ** 3 / 2 # constant for tpa: /2 because gamma = k0 * n2 + tpa / 2
heig = hevb / Eig  * 1e12         # normalisation for tpa function
x = hev * lux / lam / Eig     # 'x' for F(x) @ 1.55um
F = 0.25 * (1 / x ** 2 - 1 / x ** 3 + 0.25 / x ** 4)
F1 = 0.25 * (-2 / x ** 3 + 3 / x ** 4 - 1 / x ** 5)  # derivatives of F(x) @ 1.55um for Taylor series
F2 = 0.25 * (6 / x ** 4 - 12 / x ** 5 + 5 / x ** 6)
F3 = 0.25 * (-24 / x ** 5 + 60 / x ** 6 - 30 / x ** 7)
F4 = 0.25 * (96 / x ** 6 - 360 / x ** 7 + 210 / x ** 8) 

def tpadis(v, Rs, g, fR, gam, wm, dt, om0, h, L, dn2, dAf, dnt, tpa, Cj, dts):
    """This function is an attempt at using a Taylor expansion of tpa in terms of frequency around 1.55um"""
    # initialise arrays (26, all the same)
    b = shape(v)[0]    # length of pulses array (i.e. number of modes)
    c = shape(v)[1]    # height of array (number of time steps)
    v2 = zeros((b,c), dtype = complex)
    vf, vf2, N1, N2, N, f, A = copy(v2), copy(v2), copy(v2), copy(v2), copy(v2), copy(v2), copy(v2)
    fAb, fA, dAbdt, dAdt, iS1, iS2 = copy(v2), copy(v2), copy(v2), copy(v2), copy(v2), copy(v2)
    itpa1, itpa2, vf, d2Abdt, d2Adt = copy(v2), copy(v2), copy(v2), copy(v2), copy(v2)
    i2tpa1, i2tpa2, d3Abdt, d3Adt, d4Abdt, d4Adt = copy(v2), copy(v2), copy(v2), copy(v2), copy(v2), copy(v2)
    i3tpa1, i3tpa2, i4tpa1, i4tpa2 = copy(v2), copy(v2), copy(v2), copy(v2)

    combo = modecombo(v, Rs)
    comb2 = modecombo2(v, Rs)
    comb3 = modecombo3(v)

    # initial estimate for A(z+h,T)
    for a in range(b):
        fAb[a] = fft(ifftshift(combo[a]))
        fA[a] = fft(ifftshift(comb2[a]))
        dAbdt[a] = ifftshift(ifft(fftshift(-1j * wm) * fAb[a]))
        dAdt[a] = ifftshift(ifft(fftshift(-1j * wm) * fA[a]))
        d2Abdt[a] = ifftshift(ifft(fftshift((-1j * wm)) ** 2 * fAb[a]))
        d2Adt[a] = ifftshift(ifft(fftshift((-1j * wm)) ** 2 * fA[a]))
        d3Abdt[a] = ifftshift(ifft(fftshift((-1j * wm)) ** 3 * fAb[a]))
        d3Adt[a] = ifftshift(ifft(fftshift((-1j * wm)) ** 3 * fA[a]))
        d4Abdt[a] = ifftshift(ifft(fftshift((-1j * wm)) ** 4 * fAb[a]))
        d4Adt[a] = ifftshift(ifft(fftshift((-1j * wm)) ** 4 * fA[a]))
        iS1[a] = 1j * (1 / om0 + dn2 - dAf) * dAbdt[a]  # self-steepening term
        iS2[a] = 1j * (1 / om0 + dn2 - dAf) * conj(comb3[a]) * dAdt[a]
        itpa1[a] = con * (1j * heig * F1 * dAbdt[a])          # dispersion of tpa (F1)
        itpa2[a] = con * (1j * heig * F1 * conj(comb3[a]) * dAdt[a])
        i2tpa1[a] = con * ((1j * heig) ** 2 * F2 / 2 * d2Abdt[a])          # dispersion of tpa (F2)
        i2tpa2[a] = con * ((1j * heig) ** 2 * F2 / 2 * conj(comb3[a]) * d2Adt[a])
        i3tpa1[a] = con * ((1j * heig) ** 3 * F3 / 6 * d3Abdt[a])          # dispersion of tpa (F3)
        i3tpa2[a] = con * ((1j * heig) ** 3 * F3 / 6 * conj(comb3[a]) * d3Adt[a])
        i4tpa1[a] = con * ((1j * heig) ** 4 * F4 / 24 * d4Abdt[a])          # dispersion of tpa (F4)
        i4tpa2[a] = con * ((1j * heig) ** 4 * F4 / 24 * conj(comb3[a]) * d4Adt[a])
        sigmaf = fca(Cj, v, Rs, dts)
    
        f[a] = ifft(combo[a])  # Raman convolution
        A[a] = dt * ifftshift(fft(f[a] * g))
    
        N1[a] = 1j * real(gam) * ((1 - fR) * combo[a] + fR * A[a] + iS1[a] + iS2[a]) - imag(gam) * (combo[a] \
            + itpa1[a] + itpa2[a] + i2tpa1[a] + i2tpa2[a] + i3tpa1[a] + i3tpa2[a] + i4tpa1[a] + i4tpa2[a]) - sigmaf / 2  # Nonlinear operator
        v2[a] = v[a] * exp(h * N1[a])  # Nonlinear operator applied
        vf2[a] = (fft(ifftshift(v2[a]))) * exp(h/2 * L[a])  # Second linear operation
        v2[a] = ifftshift(ifft((vf2[a])))  # Estimate of A(z+h)

    combo = modecombo(v2, Rs)
    comb2 = modecombo2(v2, Rs)
    comb3 = modecombo3(v2)
    
    for a in range(b):
        fAb[a] = fft(ifftshift(combo[a]))
        fA[a] = fft(ifftshift(comb2[a]))
        dAbdt[a] = ifftshift(ifft(fftshift(-1j * wm) * fAb[a]))
        dAdt[a] = ifftshift(ifft(fftshift(-1j * wm) * fA[a]))
        d2Abdt[a] = ifftshift(ifft(fftshift((-1j * wm)) ** 2 * fAb[a]))
        d2Adt[a] = ifftshift(ifft(fftshift((-1j * wm)) ** 2 * fA[a]))
        d3Abdt[a] = ifftshift(ifft(fftshift((-1j * wm)) ** 3 * fAb[a]))
        d3Adt[a] = ifftshift(ifft(fftshift((-1j * wm)) ** 3 * fA[a]))
        d4Abdt[a] = ifftshift(ifft(fftshift((-1j * wm)) ** 4 * fAb[a]))
        d4Adt[a] = ifftshift(ifft(fftshift((-1j * wm)) ** 4 * fA[a]))
        iS1[a] = 1j * (1 / om0 + dn2 - dAf) * dAbdt[a]  # self-steepening term
        iS2[a] = 1j * (1 / om0 + dn2 - dAf) * conj(comb3[a]) * dAdt[a]
        itpa1[a] = con * (1j * heig * F1 * dAbdt[a])          # dispersion of tpa (F1)
        itpa2[a] = con * (1j * heig * F1 * conj(comb3[a]) * dAdt[a])
        i2tpa1[a] = con *  ((1j * heig) ** 2 * F2 / 2 * d2Abdt[a])          # dispersion of tpa (F2)
        i2tpa2[a] = con * ((1j * heig) ** 2 * F2 / 2 * conj(comb3[a]) * d2Adt[a])
        i3tpa1[a] = con * ((1j * heig) ** 3 * F3 / 6 * d3Abdt[a])          # dispersion of tpa (F3)
        i3tpa2[a] = con * ((1j * heig) ** 3 * F3 / 6 * conj(comb3[a]) * d3Adt[a])
        i4tpa1[a] = con * ((1j * heig) ** 4 * F4 / 24 * d4Abdt[a])          # dispersion of tpa (F3)
        i4tpa2[a] = con * ((1j * heig) ** 4 * F4 / 24 * conj(comb3[a]) * d4Adt[a])
        sigmaf = fca(Cj, v2, Rs, dts)
    
        f[a] = ifft(combo[a])  # Raman convolution
        A[a] = dt * ifftshift(fft(f[a] * g))
    
        N2[a] = 1j * real(gam) * ((1 - fR) * combo[a] + fR * A[a] + iS1[a] + iS2[a]) - imag(gam) * (combo[a] \
            + itpa1[a] + itpa2[a] + i2tpa1[a] + i2tpa2[a] + i3tpa1[a] + i3tpa2[a] + i4tpa1[a] + i4tpa2[a]) - sigmaf / 2  # Nonlinear operator
        N[a] = h / 2 * (N1[a] + N2[a])  # Full second-order nonlinear operator
        v[a] = v[a] * exp(N[a])  # NL operator applied
        vf[a] = fft(ifftshift(v[a])) * exp(h * L[a])  # Final linear operation
        v[a] = ifftshift(ifft(vf[a]))  # End of propagation step

    return v

def third(v, Rs, g, fR, gam, wm, dt, om0, h, L, dn2, dAf, dnt, tpa, Cj, dts, t3pa, d3pa, Cj3, nv):
    # initialise arrays
    b = shape(v)[0]    # length of pulses array (i.e. number of modes)
    c = shape(v)[1]    # height of array (number of time steps)
    v2 = zeros((b,c), dtype = complex)
    vf, vf2, N, f, A, vf, itpa1, itpa2, fAb, fA, dAbdt, dAdt, iS1, iS2, N1, N2, fA3pa, dAdt3pa, i3pa = \
        copy(v2), copy(v2), copy(v2), copy(v2), copy(v2), copy(v2), copy(v2), copy(v2), \
        copy(v2), copy(v2), copy(v2), copy(v2), copy(v2), copy(v2), copy(v2), copy(v2), copy(v2), copy(v2), copy(v2)
    qnz = 0#randn(nv) * 1e-9

    combo = modecombo(v, Rs)
    comb2 = modecombo2(v, Rs)
    comb3 = modecombo3(v)
    com3pa = modecombo3pa(v, Rs)

    # initial estimate for A(z+h,T)
    for a in range(b):
        fAb[a] = fft(ifftshift(combo[a]))
        fA[a] = fft(ifftshift(comb2[a]))
        fA3pa[a] = fft(ifftshift(com3pa[a]))
        dAbdt[a] = ifftshift(ifft(fftshift(-1j * wm) * fAb[a]))
        dAdt[a] = ifftshift(ifft(fftshift(-1j * wm) * fA[a]))
        dAdt3pa[a] = ifftshift(ifft(fftshift(-1j * wm) * fA3pa[a]))
        iS1[a] = 1j * (1 / om0 + dn2 - dAf) * dAbdt[a]  # self-steepening term
        iS2[a] = 1j * (1 / om0 + dn2 - dAf) * conj(comb3[a]) * dAdt[a]
        itpa1[a] = (1j * dnt * dAbdt[a])    # dispersion of tpa
        itpa2[a] = (1j * dnt * conj(comb3[a]) * dAdt[a])
        i3pa[a] = 1j * (d3pa * dAdt3pa[a])
        sigmaf = fca(Cj, v, Rs, dts)
        sigmaf3 = fca3(Cj3, v, Rs, dts)
    
        f[a] = ifft(combo[a])  # Raman convolution
        A[a] = dt * ifftshift(fft(f[a] * g))
    
        N1[a] = 1j * real(gam) * ((1 - fR) * combo[a] + fR * A[a] + iS1[a] + iS2[a]) - imag(gam) \
            * (combo[a] + itpa1[a] + itpa2[a]) - t3pa * com3pa[a] - i3pa[a] - sigmaf / 2 - sigmaf3 / 2 # Nonlinear operator
        v2[a] = v[a] * exp(h * N1[a])  # Nonlinear operator applied
        vf2[a] = (fft(ifftshift(v2[a]))) * exp(h/2 * L[a]) + qnz # Second linear operation
        v2[a] = ifftshift(ifft((vf2[a])))  # Estimate of A(z+h)

    combo = modecombo(v2, Rs)
    comb2 = modecombo2(v2, Rs)
    comb3 = modecombo3(v2)
    com3pa = modecombo3pa(v2, Rs)
    
    for a in range(b):
        fAb[a] = fft(ifftshift(combo[a]))
        fA[a] = fft(ifftshift(comb2[a]))
        fA3pa[a] = fft(ifftshift(com3pa[a]))
        dAbdt[a] = ifftshift(ifft(fftshift(-1j * wm) * fAb[a]))
        dAdt[a] = ifftshift(ifft(fftshift(-1j * wm) * fA[a]))
        dAdt3pa[a] = ifftshift(ifft(fftshift(-1j * wm) * fA3pa[a]))
        iS1[a] = 1j * (1 / om0 + dn2 - dAf) * dAbdt[a]  # self-steepening term
        iS2[a] = 1j * (1 / om0 + dn2 - dAf) * conj(comb3[a]) * dAdt[a]
        itpa1[a] = (1j * dnt * dAbdt[a])    # dispersion of tpa
        itpa2[a] = (1j * dnt * conj(comb3[a]) * dAdt[a])
        i3pa[a] = 1j * (d3pa * dAdt3pa[a])
        sigmaf = fca(Cj, v2, Rs, dts)
        sigmaf3 = fca3(Cj3, v2, Rs, dts)
    
        f[a] = ifft(combo[a])  # Raman convolution
        A[a] = dt * ifftshift(fft(f[a] * g))
    
        N2[a] = 1j * real(gam) * ((1 - fR) * combo[a] + fR * A[a] + iS1[a] + iS2[a]) - imag(gam) \
            * (combo[a] + itpa1[a] + itpa2[a]) - t3pa / 2 * com3pa[a] - sigmaf / 2 - sigmaf3 / 2 # Nonlinear operator
        N[a] = h / 2 * (N1[a] + N2[a])  # Full second-order nonlinear operator
        v[a] = v[a] * exp(N[a])  # NL operator applied
        vf[a] = fft(ifftshift(v[a])) * exp(h * L[a]) + qnz # Final linear operation
        v[a] = ifftshift(ifft(vf[a]))  # End of propagation step

    return v

"""Third-order method, using FT2 operator"""
def SSFMTR2(v, Rs, TR, gam, wm, dt, om0, h, L, dn2, dAf, dtpa, tpa, Cj, dts, t3pa, d3pa, d3pa2, d3pa3, Cj3, nv, nom, p0, sol, M, fine):
    # initialise arrays
    b = shape(v)[0]    # length of pulses array (i.e. number of modes)
    c = shape(v)[1]    # height of array (number of time steps)
    v2 = zeros((b,c), dtype = complex)
    vf, vf2, N, vf, itpa1, itpa2, fAb, fA, dAbdt, dAdt, iSconj, iS2, N1, N2, fAconj, dAdtconj, vin, fA3pa, dAdt3pa \
        , i3pa, dAdt3pa2, dAdt3pa3 = copy(v2), copy(v2), copy(v2), copy(v2), copy(v2), copy(v2), copy(v2), copy(v2) \
        , copy(v2), copy(v2), copy(v2), copy(v2), copy(v2), copy(v2), copy(v2), copy(v2), copy(v2), copy(v2), copy(v2) \
        , copy(v2), copy(v2), copy(v2)
        
    vin[:] = v[:]
    combo = modecombo(v, Rs)
    comb2 = modecombo2(v, Rs)
    comb3 = modecombo3(v)
    com3pa = modecombo3pa(v, Rs)
    combconj = modecomboconj(v, Rs)

    # initial estimate for A(z+h,T)
    for a in range(b):
        fAb[a] = fft(ifftshift(combo[a]))
        fA[a] = fft(ifftshift(comb2[a]))
        fA3pa[a] = fft(ifftshift(com3pa[a]))
        fAconj[a] = fft(ifftshift(combconj[a]))
        dAbdt[a] = ifftshift(ifft(fftshift(-1j * wm) * fAb[a]))
        dAdt[a] = ifftshift(ifft(fftshift(-1j * wm) * fA[a]))
        dAdt3pa[a] = ifftshift(ifft(fftshift(-1j * wm) * fA3pa[a]))
        dAdt3pa2[a] = ifftshift(ifft(fftshift(-1 * wm ** 2) * fA3pa[a]))
        dAdt3pa3[a] = ifftshift(ifft(fftshift(1j * wm ** 3) * fA3pa[a]))
        dAdtconj[a] = ifftshift(ifft(fftshift(-1j * wm) * fAconj[a]))
        iS2[a] = 1j * (2 / om0 + 2 * dn2 - 2 * dAf + 1j * TR) * conj(comb3[a]) * dAdt[a]
        iSconj[a] = 1j * (1 / om0 + dn2 - dAf + 1j * TR) * comb3[a] * dAdtconj[a]
        itpa1[a] = (1j * dtpa * dAbdt[a])    # dispersion of tpa
        itpa2[a] = (1j * dtpa * conj(comb3[a]) * dAdt[a])
        i3pa[a] = 1j * (d3pa * dAdt3pa[a] + 1j * d3pa2 * dAdt3pa2[a] - d3pa3 * dAdt3pa3[a])   # dispersion of 3pa
        sigmaf = fca(Cj, v, Rs, dts)
        sigmaf3 = fca3(Cj3, v, Rs, dts)
    
        N1[a] = 1j * real(gam) * (combo[a] + iSconj[a] + iS2[a]) - imag(gam) \
            * (combo[a] + itpa1[a] + itpa2[a]) - t3pa * com3pa[a] - i3pa[a] - sigmaf / 2 \
             - sigmaf3 / 2 # Nonlinear operator
        v2[a] = v[a] * exp(h * N1[a])  # Nonlinear operator applied
        vf2[a] = (fft(ifftshift(v2[a]))) * exp(h/2 * L[a]) # Second linear operation
        v2[a] = ifftshift(ifft((vf2[a])))  # Estimate of A(z+h)

    combo = modecombo(v2, Rs)
    comb2 = modecombo2(v2, Rs)
    comb3 = modecombo3(v2)
    com3pa = modecombo3pa(v2, Rs)
    combconj = modecomboconj(v2, Rs)
    
    for a in range(b):
        fAb[a] = fft(ifftshift(combo[a]))
        fA[a] = fft(ifftshift(comb2[a]))
        fA3pa[a] = fft(ifftshift(com3pa[a]))
        fAconj[a] = fft(ifftshift(combconj[a]))
        dAbdt[a] = ifftshift(ifft(fftshift(-1j * wm) * fAb[a]))
        dAdt[a] = ifftshift(ifft(fftshift(-1j * wm) * fA[a]))
        dAdt3pa[a] = ifftshift(ifft(fftshift(-1j * wm) * fA3pa[a]))
        dAdt3pa2[a] = ifftshift(ifft(fftshift(-1 * wm ** 2) * fA3pa[a]))
        dAdt3pa3[a] = ifftshift(ifft(fftshift(1j * wm ** 3) * fA3pa[a]))
        dAdtconj[a] = ifftshift(ifft(fftshift(-1j * wm) * fAconj[a]))
        iS2[a] = 1j * (2 / om0 + 2 * dn2 - 2 * dAf + 1j * TR) * conj(comb3[a]) * dAdt[a]
        iSconj[a] = 1j * (1 / om0 + dn2 - dAf + 1j * TR) * comb3[a] * dAdtconj[a]
        itpa1[a] = (1j * dtpa * dAbdt[a])    # dispersion of tpa
        itpa2[a] = (1j * dtpa * conj(comb3[a]) * dAdt[a])
        i3pa[a] = 1j * (d3pa * dAdt3pa[a] + 1j * d3pa2 * dAdt3pa2[a] - d3pa3 * dAdt3pa3[a])   # dispersion of 3pa
        sigmaf = fca(Cj, v2, Rs, dts)
        sigmaf3 = fca3(Cj3, v, Rs, dts)
    
        N2[a] = 1j * real(gam) * (combo[a] + iSconj[a] + iS2[a]) - imag(gam) \
            * (combo[a] + itpa1[a] + itpa2[a]) - t3pa * com3pa[a] - i3pa[a] - sigmaf / 2 \
             - sigmaf3 / 2 # Nonlinear operator
        N[a] = h / 2 * (N1[a] + N2[a])  # Full second-order nonlinear operator
        v[a] = v[a] * exp(N[a])  # NL operator applied
        vf[a] = fft(ifftshift(v[a])) * exp(h * L[a]) # Final linear operation
        v[a] = ifftshift(ifft(vf[a]))  # End of propagation step
    
    """Error checking for soliton propagation"""
    if sol:
        rel = zeros(nv)
        cond = 1
        if p0 < 4:
            cond = p0 / 4
        for num in range(nv):
            if abs(v[0][num]) >= cond and abs(vin[0][num]) >= cond:
                rel[num] = abs(1 - abs(v[0][num]) ** 2 / abs(vin[0][num]) ** 2)
        if nom > M - 2:
            print(' Maximum relative error is', '%.2e' % max(rel))
    
        if max(rel) > 4.5e-6 and fine and nom < M / 10:
            print('Whoops! Minimum relative error between iterations is', '%.2e' % max(rel), 'at', '%.1f' % (nom / M * 100) \
                  , '%')
            fine = False
        elif max(rel) < 4.5e-6 and fine and nom > M / 10:
            print(' Soliton is fine')
            fine = False
        
    return v, fine
    
def progress(M, num):
    """Takes in no. steps in loop and gives progress in 10% increments"""
    ten = int(M / 10)
    if (num + 1) % ten == 0:
                per = int((num + 1)/ten * 10)
                sys.stdout.write ('\b\b\b%s' % per) #\b character backspace-deletes previous output
                sys.stdout.flush()
                
def modecombo(vs, Rs): 
    """takes in array of pulses ('v's), one for each mode, and corresponding overlaps ('R's) where 1/R11 = Aeff11 etc"""
    b = shape(vs)[0]    # length of pulses array (i.e. number of modes)
    c = shape(vs)[1]    # height of array (number of time steps)
    vso = zeros((b, c))    # set output pulses array to shape of input
    for a in range(b):
        vso[a] = Rs[a][a] * abs(vs[a]) ** 2 # self-term in nonlinear operator
        for d in range(b):
            if d != a:
                vso[a] += 2 * Rs[a][d] * abs(vs[d]) ** 2    # cross-terms in nonlinear operator
    return vso

def modecombo3pa(vs, Rs):
    """Self and cross terms in generated by 3pa (i.e. proportional to mod^4)"""
    b = shape(vs)[0]    # length of pulses array (i.e. number of modes)
    c = shape(vs)[1]    # height of array (number of time steps)
    vso = zeros((b, c))    # set output pulses array to shape of input
    for a in range(b):
        vso[a] = Rs[a][a] ** 2 * abs(vs[a]) ** 4 # self-term in nonlinear operator
        for d in range(b):
            if d != a:
                vso[a] += 2 * Rs[a][d] ** 2 * abs(vs[d]) ** 4    # cross-terms in nonlinear operator
    return vso
    
def modecombo4pa(vs, Rs):
    """Self and cross terms in generated by 3pa (i.e. proportional to mod^4)"""
    b = shape(vs)[0]    # length of pulses array (i.e. number of modes)
    c = shape(vs)[1]    # height of array (number of time steps)
    vso = zeros((b, c))    # set output pulses array to shape of input
    for a in range(b):
        vso[a] = Rs[a][a] ** 3 * abs(vs[a]) ** 6 # self-term in nonlinear operator
        for d in range(b):
            if d != a:
                vso[a] += 2 * Rs[a][d] ** 3 * abs(vs[d]) ** 6    # cross-terms in nonlinear operator
    return vso

def modecombo2(vs, Rs):
    """As modecombo but just pulse amplitude (not mod or squared) to use in conjugate operator"""
    b = shape(vs)[0]    # length of pulses array (i.e. number of modes)
    c = shape(vs)[1]    # height of array (number of time steps)
    vso = zeros((b, c), dtype = complex)    # set output pulses array to shape of input
    for a in range(b):
        vso[a] = vs[a] * Rs[a][a] # self-term in nonlinear operator
        for d in range(b):
            if d != a:
                vso[a] += 2 * vs[d] * Rs[a][d]   # cross-terms in nonlinear operator
    return vso

def modecombo3(vs):
    """As modecombo2 but without overlap integrals"""
    b = shape(vs)[0]    # length of pulses array (i.e. number of modes)
    c = shape(vs)[1]    # height of array (number of time steps)
    vso = zeros((b, c), dtype = complex)    # set output pulses array to shape of input
    for a in range(b):
        vso[a] = vs[a] # self-term in nonlinear operator
        for d in range(b):
            if d != a:
                vso[a] += 2 * vs[d]   # cross-terms in nonlinear operator
    return vso

def modecomboconj(vs, Rs): # takes in array of pulses ('v's), one for each mode, and corresponding overlaps ('R's) where 1/R11 = Aeff11 etc
    b = shape(vs)[0]    # length of pulses array (i.e. number of modes)
    c = shape(vs)[1]    # height of array (number of time steps)
    vso = zeros((b,c), dtype = complex)    # set output pulses array to shape of input
    for a in range(b):
        vso[a] = conj(vs[a]) * Rs[a][a] # self-term in nonlinear operator
        for d in range(b):
            if d != a:
                vso[a] += 2 * conj(vs[d]) * Rs[a][d]   # cross-terms in nonlinear operator
    return vso

def fca(Cj, vs, Rs, dt):    # Cj is constants and also contains sigma * (1 + 1j * muf)
    """Calculates free carriers from cumulative integral of pulse in time"""
    b = shape(vs)[0]    # length of pulses array (i.e. number of modes)
    c = shape(vs)[1]    # height of array (number of time steps)
    Nc = zeros(c, dtype = complex)    # set output pulses array to shape of input
    for u in range(b):
        Nc += Cj * Rs[u][u] ** 2 * cumtrapz(abs(vs[u]) ** 4, None, dt, -1, 0)
        for d in range(b):
            if d != u:
                Nc += 2 * Cj * Rs[u][u] * Rs[u][d] * cumtrapz(abs(vs[u] * vs[d]) ** 2, None, dt, -1, 0)
    return Nc
    
def fca3(Cj, vs, Rs, dt):    # Cj is constants and also contains sigma * (1 + 1j * muf)
    """Calculates free carriers from cumulative integral of pulse in time"""
    b = shape(vs)[0]    # length of pulses array (i.e. number of modes)
    c = shape(vs)[1]    # height of array (number of time steps)
    Nc = zeros(c, dtype = complex)    # set output pulses array to shape of input
    for u in range(b):
        Nc += Cj * Rs[u][u] ** 3 * cumtrapz(abs(vs[u]) ** 6, None, dt, -1, 0)
        for d in range(b):
            if d != u:
                Nc += 2 * Cj * (Rs[u][u] * Rs[u][d]) ** 1.5 * cumtrapz(abs(vs[u] * vs[d]) ** 3, None, dt, -1, 0)
    return Nc    

def fca4(Cj, vs, Rs, dt):    # Cj is constants and also contains sigma * (1 + 1j * muf)
    """Calculates free carriers from cumulative integral of pulse in time"""
    b = shape(vs)[0]    # length of pulses array (i.e. number of modes)
    c = shape(vs)[1]    # height of array (number of time steps)
    Nc = zeros(c, dtype = complex)    # set output pulses array to shape of input
    for u in range(b):
        Nc += Cj * Rs[u][u] ** 4 * cumtrapz(abs(vs[u]) ** 8, None, dt, -1, 0)
        for d in range(b):
            if d != u:
                Nc += 2 * Cj * (Rs[u][u] * Rs[u][d]) ** 2 * cumtrapz(abs(vs[u] * vs[d]) ** 4, None, dt, -1, 0)
    return Nc

def sech(x):
    """Defines sech function as 1/cosh. This is then vectorized below to take in arrays"""
    if -700 < x < 700:
        return cosh(x) ** (-1)
    else:
        return 0
    
sech = vectorize(sech, otypes = [float])