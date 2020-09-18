import logging
from pymcr.constraints import *
from pymcr.regressors import OLS, NNLS
from pymcr.mcr import McrAR
from funcs import *
import sys
import matplotlib.pyplot as plt
from scipy.signal import resample
from math import floor

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
chemical_rank_tolerance = float(dct['chemical rank tolerance'])
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
    print('resampling')
    A = resample(A, num, axis=0)
    time_per_scan /= num / num0
    print('upsampled from {0} to {1} scans'.format(number_scans, A.shape[0]))

#import sys; sys.exit()
A = A.transpose()
print('period averaging')
B = period_average(A,
                   time_per_period=time_per_period,
                   time_per_scan=time_per_scan)
scans_per_period = floor(time_per_period / time_per_scan)
n_periods = round(A.shape[1]/scans_per_period)
print('period-averaged over {} periods'.format(n_periods))
print('from {0} scans to {1} scans'.format(A.shape[1], scans_per_period))

if dct['animate time dependent data'] == 'True':

    import matplotlib
    from matplotlib.animation import FuncAnimation

    matplotlib.rcParams['animation.bitrate'] = 2000
    fig, ax = plt.subplots(nrows=1, ncols=1)
    l1, = ax.plot(x, B[:, 0], label='avg')
    l2, = ax.plot(x, A[:, 0], label='transient')
    yl = np.min(np.concatenate((np.ravel(A), np.ravel(A))))
    yh = np.max(np.concatenate((np.ravel(A), np.ravel(A))))
    ax.set_ylim(yl, yh)
    ax.set_ylabel('Spectral Intensity (a.u.)')
    ax.set_xlabel('Wavelength in wavenumbers')
    ax.legend(loc='upper center')

    def init():
        l1.set_ydata([np.nan] * len(A[0]))
        l2.set_ydata([np.nan] * len(A[0]))
        return l1, l2

    def animate(i):
        l1.set_ydata(B[:, i % scans_per_period])
        l2.set_ydata(A[:, i])
        title_str = """{0}/{1} (resampled) period averaged scans
number periods = {2}
time per scan = {3:.2f}s""".format(i % scans_per_period
        , scans_per_period
        , n_periods 
        , time_per_scan)
        ax.set_title(title_str)

        return l1, l2

    interval = max((round(60*1000/scans_per_period), 50))
    ani = FuncAnimation(fig, animate,
                        frames=scans_per_period,
                        init_func=init,
                        interval=interval,
                        blit=False,
                        save_count=50)
    plt.show()
    #import sys; sys.exit()

# phase sensitive detection transform
print('calculating time to phase space transform')
C = psd_transform(B)[..., 0]
print('calculated from time space scans {} to 360 phase space'.format(scans_per_period)) 
# baseline subtraction by piecewise linear functions
pin_indices = [np.argmin(abs(x - pin)) for pin in pin_wavelengths]
print('baseline subtracting')
D = baseline_subtract(C, pin_indices=pin_indices)
print('baseline subtracted at pin wavelengths {}'.format(pin_wavelengths))
if dct['plot baseline subtracted spectra'] == 'True':
    plt.figure()
    plt.plot(D.transpose())
    plt.show()

# NMF from sklearn
#from sklearn.decomposition import non_negative_factorization
#if dct['number components'] in ['calculate', '']:
#    n_components = min(A.shape)
#else:
#    n_components = int(dct['number components'])
#
#C,ST,niter = non_negative_factorization(A,
#        n_components=n_components,
#        init='nndsvda',
#        solver='cd')
#
#import sys; sys.exit()

# singular value decomposition
print('calculating singular value decomposition')
svd, l = snglr_vl_dcmpstn(D, k_range=k_range)
print('calculated svd to number components {}'.format(int(dct['maximum number of components'])))

l = np.array(l)

if dct['number components'] in ['calculate', '']:
    ncomponents = np.argmin(l/l.max() > chemical_rank_tolerance) + 1  
    # by default argmin returns the first index of the minimum value when it appears more than once
else:
    ncomponents = int(dct['number components'])

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

#print("Approximating by {0} components".format(i))

u, s, vt = svd
ST_guess = vt[:ncomponents]

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
print('calculating MCR-ALS')
mcrar.fit(D, ST=ST_guess, verbose=True)

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
