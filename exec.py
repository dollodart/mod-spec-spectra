import logging
from pymcr.constraints import *
from pymcr.regressors import OLS, NNLS
from pymcr.mcr import McrAR
from funcs import *
import sys
import matplotlib.pyplot as plt
from scipy.signal import resample

# read user input parameters
dct = read_parameters()
number_wavelengths = dct['number wavelengths']
number_scans = dct['number scans']
dr = dct['data directory']
wavelength_range = dct['wavelength range (1/cm)']
pin_wavelengths = np.array([float(x)
                            for x in dct['pin wavelengths'].split(' ')])
xbounds = dct['wavelength bounds (1/cm)'].split(' ')
xlb, xub = float(xbounds[0]), float(xbounds[1])
k_range = range(1, int(dct['maximum number of components']) + 1)
error_tolerance = float(dct['error tolerance'])
max_iter = int(dct['maximum number of iterations'])

if number_wavelengths in [
        'calculate',
        ''] or number_scans in [
            'calculate',
            ''] or wavelength_range in [
                'calculate',
        '']:
    from pathlib import Path
    files = list(Path(dr).iterdir())
    files = [str(i) for i in files if i.suffix == '.csv']
    number_scans = len(files)
    x = np.genfromtxt(files[0], delimiter=',').T[0]
    number_wavelengths = len(x)
else:
    number_wavelengths = int(number_wavelengths)
    number_scans = int(number_scans)
    x = np.arange(int(wavelength_range[0]), int(wavelength_range[1]) + 1)

# load and reshape data
data = load(dr)
A = data.reshape((number_scans, number_wavelengths))
# A = data.reshape((number_wavelengths,number_scans))
# truncation
bl = np.logical_and(xlb < x, x < xub)
x, A = x[bl], A[:, bl]
number_wavelengths = len(x)
# data averaging
time_per_scan = float(dct['time per scan (s)'])
# if float and not int, different behavior
time_per_period = float(dct['time per period (s)'])
# resampling time axis
if dct['resample'] == 'True':
    rs_ref = float(dct['resampling refinement'])
    from math import ceil
    num0 = A.shape[0]
    # uses the Fourier Transform algorithm which is most efficient when the
    # number of samples is a power of 2
    num = 2**ceil(np.log2(rs_ref * number_wavelengths))
    # this efficiency improvement is so pronounced that increasing the number
    # of data points to the next power of 2 from a desired number is always
    # (in practice) more efficient
    A = resample(A, num, axis=0)
    time_per_scan /= num / num0

A = A.transpose()

if dct['animate time dependent data'] == 'True':
    full_data = A.copy()
    A = period_average(A,
                       time_per_period=time_per_period,
                       time_per_scan=time_per_scan)
    import matplotlib
    from matplotlib.animation import FuncAnimation

    from math import floor
    scans_per_period = floor(time_per_period / time_per_scan)

    matplotlib.rcParams['animation.bitrate'] = 2000
    fig, ax = plt.subplots(nrows=1, ncols=1)
    l1, = ax.plot(x, A[:, 0], label='avg')
    l2, = ax.plot(x, full_data[:, 0], label='transient')
    yl = np.min(np.concatenate((np.ravel(full_data), np.ravel(A))))
    yh = np.max(np.concatenate((np.ravel(full_data), np.ravel(A))))
    ax.set_ylim(yl, yh)
    ax.legend(loc='upper center')

    def init():
        l1.set_ydata([np.nan] * len(A[0]))
        l2.set_ydata([np.nan] * len(full_data[0]))
        return l1, l2

    def animate(i):
        l1.set_ydata(A[:, i % scans_per_period])
        l2.set_ydata(full_data[:, i])
        ax.set_title("{0}/{1} and {2}/{3}".format(i %
                                                  scans_per_period, scans_per_period, i, number_scans))
        return l1, l2

    ani = FuncAnimation(fig, animate,
                        frames=number_scans,
                        init_func=init,
                        interval=200,
                        blit=False,
                        save_count=50)
    plt.show()

else:
    A = period_average(A,
                       time_per_period=time_per_period,
                       time_per_scan=time_per_scan)

# phase sensitive detection transform
A = psd_transform(A)[..., 0]
# baseline subtraction by piecewise linear functions
pin_indices = [np.argmin(abs(x - pin)) for pin in pin_wavelengths]
A = baseline_subtract(A, pin_indices=pin_indices)
if dct['plot baseline subtracted spectra'] == 'True':
    plt.figure()
    plt.plot(A.transpose())
    plt.show()
# singular value decomposition
svd, l = snglr_vl_dcmpstn(A, k_range=k_range)
if dct['plot SVD vectors'] == 'True':
    fig, axs = plt.subplots(nrows=1, ncols=3)
    axs[0].plot(svd[0])
    axs[1].plot(svd[1])
    axs[2].plot(svd[2].T)
    plt.show()
if dct['plot SVD error'] == 'True':
    plt.figure()
    plt.xlabel('Number SVD Components')
    plt.ylabel('Error in SVD Approximation')
    plt.semilogy(k_range, l)
    plt.show()

if dct['number components'] in ['calculate', '']:
    i = np.argmin(np.array(l) > error_tolerance) + 1  
    # by default argmin returns the first index of the minimum value when it appears more than once
else:
    i = int(dct['number components'])

#print("Approximating by {0} components".format(i))

u, s, vt = svd
ST_guess = vt[:i]

# multivariate curve resolution
c_constraints = []
st_constraints = []

# likely a more elegant way exists
c_l = dct['concentration constraints'].split(' ')
st_l = dct['spectrum constraints'].split(' ')
if 'ConstraintNonneg' in c_l:
    c_constraints.append(ConstraintNonneg())
if 'ConstraintCumsumNonneg' in c_l:
    c_constraints.append(ConstraintCumsumNonneg())
if 'ConstraintNorm' in c_l:
    c_constraints.append(ConstraintNorm())
if 'ConstraintNonneg' in st_l:
    st_constraints.append(ConstraintNonneg())
if 'ConstraintCumsumNonneg' in st_l:
    st_constraints.append(ConstraintCumsumNonneg())
if 'ConstraintNorm' in st_l:
    st_constraints.append(ConstraintNorm())
# print(c_constraints,st_constraints)
logger = logging.getLogger('pymcr')
logger.setLevel(logging.DEBUG)
stdout_handler = logging.StreamHandler(stream=sys.stdout)
stdout_format = logging.Formatter('%(message)s')
stdout_handler.setFormatter(stdout_format)
logger.addHandler(stdout_handler)

mcrar = McrAR(
    max_iter=max_iter,
    st_regr=NNLS(),
    c_regr=OLS(),
    c_constraints=[],
    st_constraints=st_constraints)
mcrar.fit(A, ST=ST_guess, verbose=True)

np.savetxt('out/wavelengths.txt', x)
np.savetxt('out/C-opt.txt', mcrar.C_opt_)
if dct['smooth spectra'] == 'True':
    args = (int(dct['smoothing window length']), str(dct['smoothing window']))
    S = np.apply_along_axis(smooth, 1, mcrar.ST_opt_, *args).T
    np.savetxt('out/S-opt.txt', S)
else:
    S = mcrar.ST_opt_.T
    np.savetxt('out/S-opt.txt', S)

if dct['plot MCR ALS results'] == 'True':
    plt.figure(figsize=(6, 4))
    plt.subplot(211)
    plt.plot(mcrar.C_opt_)
    plt.xlabel('Phase angle in deg.')
    plt.ylabel('Concentration (arb. units)')
    plt.subplot(212)
    plt.xlabel('Wavelength in wavenumbers')
    plt.ylabel('Spectral Intensity (arb. units)')
    plt.plot(x, S)
    plt.tight_layout()
    plt.show()
