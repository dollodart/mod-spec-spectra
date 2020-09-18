# Overview

This script calculates principal components from transient spectra taken during modulation excitation spectroscopy experiments. The idea is that different chemical species will have, in addition to different component spectra and concentrations, different transient responses to modulation of reactant concentrations. By transforming to the phase domain and finding the distinct species spectra and concentration by a method such as alternating least squares based off their different transient and therefore phase behavior, insight into the reaction mechanisms can be obtained. An understanding of the relative concentrations and intrinsic spectral intensities of pure component species can assist in the principal component determination.

This program was developed for users not experienced with Python and working in a windows environment and installation and usage instructions reflect that.

# Installation

This script is written in Python 3.7. Make sure your python version is at least 3.7. Since Python 2 was held for so long that sometimes even if version 3 is installed, by default `python` is python 2.7 and the version 3 needs to be specified as `python3`. To check the version, in the IPython shell command line type `python -V` and `python3 -V`. This script has several python packages it uses, which are called dependencies. They are 

- sys
- pickle
- math
- pathlib
- logging

- scipy
- numpy
- matplotlib
- natsorted
- pymcr (stylized pyMCR)

Each dependency has a version or range of versions required. I developed the script with the most recent version of these packages as of 2019-10-24. There are many online tutorials for downloading and installing python packages, here is one:

https://www.earthdatascience.org/courses/earth-analytics-python/python-open-science-toolbox/install-and-import-python-packages/

If the package source is irrelevant and the package can be installed in the system python directory, type the command `conda install <package name>` where you substitute the package name for the angle bracketed term. Note, you may have to do something like `conda3 install <package name>` in order to have the python 3 package version download.

# Usage

Assuming you have the program and data and know their respective paths

1. Use wordpad to change the parameters as needed for your data in the 'user-inputs.txt' file in the program root directory.
2. Open up the Ipython shell. This can usually be done by typing IPython in the Start search bar and clicking the shortcut.
3. In the Ipython shell, navigate to the program root directory.
4. In the command line of the Ipython shell, type `python3 exec.py`. This calls the python interpreter to execute the script.
5. Depending on your parameter choice some interactive pop-ups may appear. Close the pop-ups in order to let the script proceed to completion.
6. The final data are stored in the 'out/' directory in the program root directory as csv files.

# Description of parameters

- root directory: This is the program root directory, wherever you have the script and the associated directories. Note that Windows accepts paths with directories delimited by / (POSIX standard) as well as \ (DOS convention). It must have a trailing slash to indicate a directory, rather than a regular file. For example, if the program root directory is C:\Programs\spectra, you would input C:/Programs/spectra/.
- data directory: This is the directory which contains the csv files. The same conditions apply as those on the root directory.
- number scans: This is the total number of scans. When given 'calculate' this is calculated as the number of files in the data directory.
- number wavelengths: This is the total number of wavelengths at which the spectrum was taken. When given 'calculate' this is calculated from the data provided.
- wavelength range (1/cm): This is the range of wavelengths (in Python, a range is an iterable, like list(range(1,10))=[1,2,3,4,5,6,7,8,9]). When given 'calculate' this is taken from the data provided.
- maximum number of components: This the maximum number of components the SVD decomposition will calculate.
- number components: This is the ideal number of components. When given 'calculate' it chooses the minumum number of components such that the error tolerance is satisfied.
- time per period (s): This is the modulation period used in seconds.
- time per scan (s): This is the time between acquisitions in seconds.
- pin wavelengths: This the list of wavelengths at which to interpolate piecewise linear functions for the baseline subtraction
- wavelength bounds (1/cm): This is the range of wavelengths desired. The data outside of these bounds is omitted.
- error tolerance: This is the maximum allowed error in the data matrix decomposition for the SVD guess. It is related to the experimental noise expected. It has the same arbitrary units of intensity as the data matrix.
- resample: This says whether or not to resample the data by a Fourier transform method.
- resampling refinement: This is the minimum factor of the refinement of the upsampling done. Generally the refinement is higher since the fast fourier transform requires the number of samples to be a power of 2, and the next highest power of 2 from that calculated by the refinement is used.
- smooth spectra: This says whether or not to apply smoothing to the spectra after the MCR-ALS.
- smoothing window: This is the type of smoothing window, which are how points in the window are weighed in the smoothing (usually points in the center are weighed more than those on the edges, but in the case of flat there is equal weighting and therefore greater smoothing). See references for details. Allowed values are 'flat','hanning','hamming','bartlett','blackman'.
- smoothing window length: The number of points included in the smoothing window. This must be smaller than the number of data points.
- concentration constraints: Constraints to be applied to the concentration matrix in the MCR-ALS iteration. Can be values of 'ConstraintNonneg', 'ConstraintCumsumNonneg', or 'ConstraintNorm'
- spectrum constraints: Constraints to be applied to the spectrum matrix in the MCR-ALS iteration. Can be values of 'ConstraintNonneg', 'ConstraintCumsumNonneg', or 'ConstraintNorm'
- maximum number of iterations: maximum number of iterations allowed for the MCR-ALS algorithm
- animate time dependent data: Animates the raw data and the period averaged data to see if there is agreement between the two
- plot SVD error: Plots the error of the SVD as a function of number of components
- plot SVD vectors: Plots the component spectra and concentrations based on the
- plot baseline subtracted spectra: Plots the PES-MSD transformed and baseline subtracted spectra as "isophases"
- plot MCR ALS results: Plots the final results, that is, concentrations and spectra of deconvoluted components

# Notes on types for user inputs
- Booleans will evaluate True for 'True' and False otherwise
- For values to be automatically calculated, either assign no value '', or the value 'calculate'
- Ranges are supplied as two numbers separated by a space, e.g., '10 20'
- Lists are supplied with spaces between elements, e.g., '10 20 30 40 50'
- Excess whitespace in user inputs file is stripped from the left and right, e.g., ' 10' or '10 ' will turn to '10'

# Example user input file
A set of example values for each parameter is given below (it also gives the type of data, special cases which can be given the 'calculate' value are denoted with a pipe, in an actual input file only one of the inputs on either side of the pipe can be input):

root directory,/path/to/my/root/dir/
data directory,/path/to/my/data/dir/
number scans,400|calculate
number wavelengths,300|calculate
wavelength range (1/cm),100 200|calculate
maximum number of components,10|calculate
number components,5|calculate
time per period (s),300.
time per scan (s),5.0
pin wavelengths,100 200 300 400
wavelength bounds (1/cm),100 300
error tolerance,0.01
resample,True
resampling refinement,2.0
smooth spectra,True
smoothing window,flat
smoothing window length,20
concentration constraints,ConstraintNorm
spectrum constraints,ConstraintCumsumNonneg
maximum number of iterations,100
animate time dependent data,False
plot SVD error,False
plot SVD vectors,False
plot baseline subtracted spectra,False
plot MCR ALS results,True

# Description of Program

## Acronyms

- mes=modulation excitation spectroscopy
- psd=phase sensitive detection
- mcr=multivariate curve resolution
- als=alternating least squares
- lof=local outlier factor
- svd=singular value decomposition
- pca=principal components analysis

## Definitions 

Chemical Rank: The rank of a data matrix found by the
number of singular values which are greater than some cut-off value
determined by experimental error.

## Data Transformation

Begin with time resolved spectra in a matrix with wavelength as rows and
time as columns, intensity as entries. It may be necessary, if there
are multiple data sets which are to be processed at once by augmented
matrices, to do interpolation to obtain equal values of the tested
conditions. For the following it is assumed there is no augmented matrix
and only one data matrix is being used.

1. For the experimental data used, to eliminate noise, several periods
of data are taken, and the signal at a given time modulo the period is
averaged to reduce the experimental noise. For example, if the period
is 2 seconds, and 10 seconds and 5 periods are taken, then the data
point at 1.3 seconds is averaged with that at 3.3, 5.3, 7.3, and 9.3
seconds. This averaging requires resampling (effectively interpolation)
since the data is generally not collected at a frequency which is an
integer multiple of the modulation frequency.

2. Perform a phase-sensitive transform to the phase angle space
to obtain a matrix with wavelength as rows and mode number as columns,
intensity as entries.

Because some of the further operations require evaluating first and
second derivatives, it may be necessary to smooth the data by a method
such as spline or Whittaker smoothing. 

3. Evaluate the SVD decomposition of the data matrix. If the original
data matrix has features with different time variations, that will be
preserved in the transform to the phase angle space, and the SVD will
give as many distinct components as there are different time varying
features in the data matrix. Since intermediates and products may not
oscillate at the modulation frequency used for the reactants, to first
approximation this is the number of chemical species present. But the
SVD will give a spectrum for each component of arbitrary complexity
(number of features), and only discriminates based on time/phase angle
variation. The SVD decomposition is often written USV^T, and in this
case the U matrix has as column vectors are the time/phase angle varying
concentrations, the V^T has as column vectors the spectra, and the S
has as the diagonal entries the weight of that species to the observed
spectroscopic signal.

4. For each resulting spectrum in the SVD decomposition, make an
initial guess for the spectrum of each component by peak assignments
in the literature. This will be greater than the number of guessed
components in the SVD since the component spectrum can be arbitrarily
complicated, the only cause of different components being differences in 
transient response. Then initialize the spectrum matrix transpose by
stacking the component spectra. The SVD decomposition gives the chemical
rank by the number of non-zero entries in the diagonal matrix which are
above the error tolerance. Note S is determined up to a multiplicative
constant unless U and V^T are normalized which is relevant when imposing
an error tolerance. In this implementation the SVD decomposition is used
directly as the initial guess for step 5.

5. Perform a MCR-ALS fit on the convoluted data matrix using the
component spectrum matrix found in step 4 as an initial guess. This
gives a regressed concentration matrix in addition to the spectrum
matrix. The concentration matrix has modes as rows and wavelengths
as columns. This might be interpreted as each component having some
different effective concentration at higher harmonics. In this
implementation only the fundamental mode is used, and so this is a
vector of concentrations at the fundamental mode (the modulation
frequency).

Note the SVD decomposition is used in forming the psuedo inverse to
solve the general least squares problem. Therefore when the SVD is input
to a multivariate curve resolution alternating least squares as an
initial value, the convergence should be immediate since it is already a
solution to the least squares problem. It remains to be implemented to
allow the user to specify the pure component spectra based off the SVD
so that the MCR-ALS gives a different result than the SVD.

The overview talks about principal component analysis, see, e.g., https://stats.stackexchange.com/questions/134282/relationship-between-svd-and-pca-how-to-use-svd-to-perform-pca. The SVD contains the same and more information than PCA.