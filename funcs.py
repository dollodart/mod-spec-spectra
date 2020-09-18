from scipy.sparse.linalg import svds
from math import floor
from pathlib import Path
import numpy as np
from natsort.natsort import natsort
import pickle

with open('./user-inputs.txt', 'r') as read_file:
    for line in read_file.readlines():
        line = line.rstrip('\n')
        line = line.split(',')
        line = [i.strip() for i in line]
        if line[0] == 'root directory':
            ROOT_DIR = line[1]


def compress(dr):
    """

    Takes directory containing xy spectra in csv format as string and
    compresses to a binary format using the python pickle library.

    The data is saved with the wavelength index fastest and time index
    index slowest (the default option for numpy ravel is C-style, or row
    major, order).

    """

    files = list(Path(dr).iterdir())
    files = [str(i) for i in files if i.suffix == '.csv']
    files = natsort.humansorted(files)
    l = []
    for f in files:
        y = np.genfromtxt(f, delimiter=',').T[1]
        l.append(y)
    with open(ROOT_DIR + "data/{0}.pickle".format(dr.split('/')[-2]), 'wb') as write_file:
        pickle.dump(np.ravel(l), write_file)
    return None


def load(dr):
    """

    Loads data, first checking if the compressed form exists, otherwise
    parsing the data directory.

    """

    save_name = ROOT_DIR + "data/{0}.pickle".format(dr.split('/')[-2])
    if Path(save_name).exists():
        with open(save_name, 'rb') as binary_read_file:
            return pickle.load(binary_read_file)
    else:
        compress(dr)
        with open(save_name, 'rb') as binary_read_file:
            return pickle.load(binary_read_file)


def read_parameters():
    try:
        with open('user-inputs.txt', 'r') as read_file:
            lines = read_file.readlines()
            kv = []
            for counter, line in enumerate(lines):
                line = line.rstrip('\n')  # remove newline
                line = line.split(',')
                line = [i.strip() for i in line]
                if len(line) == 2:
                    kv.append(line)
                else:
                    print(
                        'omitted line {0} at line number {1}'.format(
                            ''.join(line), counter))
            return dict(kv)

    except Exception as e:
        print('Exception (likely the parameters file is not formatted correctly')
        print(e)


def period_average(data, time_per_period, time_per_scan):
    """

    Takes data matrix with rows corresponding to the independent
    spectroscopy variable (wavelength) and columns corresponding to
    times and averages spectra taken at the same point in different
    periods.

    """

    number_wavelengths, number_scans = data.shape
    scans_per_period = floor(time_per_period / time_per_scan)
    number_periods = floor(number_scans / scans_per_period)
    remainder = time_per_period / time_per_scan - time_per_period // time_per_scan
    if remainder < 0.5:
        sgn = -1
    else:
        sgn = 1
        remainder = 1. - remainder

    # for floats, x % y != x/y - x // y

    avg_data = np.zeros((number_wavelengths, scans_per_period))

    for i in range(number_periods):
        offset = round(remainder * i)
        lb = i * scans_per_period + sgn * offset
        ub = (i + 1) * scans_per_period + sgn * offset
        intrvl = range(lb, ub)
        avg_data += data[:, intrvl]

    avg_data /= number_periods
    return avg_data


def ntgrtn_cffcnts(n, order=1):
    """

    Returns integration coefficients depending on desired order of
    accuracy. Does not normalize them.

    """

    if order == 1:
        coeffs = np.ones(n)
        coeffs[0], coeffs[-1] = 0.5, 0.5
    else:
        raise Warning('Order > 1 not yet supported')
    return coeffs


def psd_transform(A, phis=np.linspace(0, 2 * np.pi, 360), ks=np.arange(1, 2)):
    """

    The provided matrix must correspond to one period.

    """
    nr, nc = A.shape
    nphis = len(phis)
    args = np.linspace(0, 1, nc) * 2 * np.pi
    args = np.outer(args, ks)
    args = np.tile(args, reps=(nphis, 1, 1)) + phis[:, np.newaxis, np.newaxis]
    coeffs = ntgrtn_cffcnts(nc) * 2. / nc
    integral = A * coeffs
    return integral @ np.sin(args)


def baseline_subtract(A, pin_indices=[0, -1]):
    """

    Subtracts baseline given by piecewise linear functions of the pin
    indices.  Assumes the independent variable is evenly spaced.  Makes
    a copy of the data. An inplace substitution can be done for memory
    concerns by saving the relevant data (at each pin) before running
    the loop.

    """

    B = A.copy()

    i0 = 0
    if pin_indices[0] == 0:
        pin_indices = pin_indices[1:]
    for i in pin_indices:
        m = (A[:, i] - A[:, i0]) / (i - i0)
        dx = np.arange(i - i0)
        l = np.outer(m, dx) + A[:, i0, np.newaxis]
        B[:, i0:i] -= l
        i0 = i
    if i != A.shape[1] - 1:
        m = (A[:, -1] - A[:, i]) / (A.shape[1] - i)
        dx = np.arange(i, A.shape[1])
        l = np.outer(m, dx) + A[:, i, np.newaxis]
        B[:, i:] -= l

    return B


def snglr_vl_dcmpstn(A, k_range=range(1, 51)):
    """

    Calculates singular value decomposition provided range of number
    of singular values and quantifies error in approximation with each
    number of singular components in the normalized Frobenius norm.

    """

    u, s, vt = svds(A, max(k_range))
    l = []
    for k in k_range:
        diag = np.diag(s[:k])
        approxA = u[:, :k]@diag@vt[:k]
        E = A - approxA
        fnorm = np.linalg.norm(E)  # normalize this
        l.append(fnorm)
    return (u, s, vt), l


def smooth(x, window_len=11, window='hanning'):
    """

    Smooth the data using a window with specified length.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.

    input:
        x: the input signal
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        y: the smoothed signal

    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.

    """

    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")

    if window_len < 3:
        return x

    if window not in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError(
            "Window must be one of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")

    dct = {
        'hanning': np.hanning,
        'hamming': np.hamming,
        'bartlett': np.bartlett,
        'blackman': np.blackman}
    # s = np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]] 
    # this concatenates x with its first window length of points (reversed) and second window length of points (in original order)
    # this is bigger than the window length, and hence in the convolution
    # there is a "traveling window" around each position
    s = x

    if window == 'flat':  # moving average
        w = np.ones(window_len, 'd')
    else:
        w = dct[window](window_len)

    # discrete convolution is given for all points, these generally differ in
    # length
    y = np.convolve(w / w.sum(), s, mode='same')
    # convolution(vector1,vector2,position) = sum(vector1*np.roll(vector2,position))
    # for positions 0,1,...,min([vector1.size,vector2.size])
    # with valid casting mode supports vectors of different sizes
    return y
