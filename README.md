# UPbplot.py

[Introduction]

This is a script for calculation and visualization tool of U-Pb age
data.  It enables us to

1. plot scattered data with error ellipses on conventional
(207Pb/235U-206Pb/238U) and Terra--Wasserburg (238U/206Pb-207Pb/206Pb)
concordia diagrams, and

2. calculate one- or two-dimensional weighted mean, concordia, and
concordia-intercept ages with errors on both concordia diagrams.


[License, Author, and Copyright]

	License: Apache License, Version 2.0
	Author: Atsushi Noda
	Copyright: Geological Survey of Japan, AIST


[Citation]

Noda, Atsushi (2017) A new tool for calculation and visualization
of U-Pb age data: UPbplot.py.  Bulletin of the Geological Survey of
Japan, 68(3), 131-140, https://doi.org/10.9795/bullgsj.68.131.


[Development environment]

The script was written in Python 3.8.2 and tested on MacOSX
(Mojave).  If you use it under Python 2, please modify lines using 
print function in the script.  Comments are added to the relevant lines.

Mandatory and optional libraries are listed in [Preparation] section.


[Disclaimer]

All scripts and software are provided "as is", without any
warranty or guarantee of its usability or fitness for any purpose.
The author disclaims all liability for direct, indirect, or
consequential damages resulting from your use of the programs.


[Support]

The author does not provide any support on the usage of Python and
related libraries.  See the documentation and help resources on their
websites if you need such help.


[Preparation]

1. When you use this script (UPbplot.py) from this source code,
   additional libraries (Numpy, SciPy, Matplotlib, pandas, and so on)
   will be required.  Install them in advance.

		Numpy: https://pypi.org/project/numpy/
		Matplotlib: https://pypi.org/project/matplotlib/
		pandas: https://pypi.org/project/pandas/
		SciPy: https://pypi.org/project/scipy/
		Shapely: https://pypi.org/project/Shapely/

2. If you prefer to using Qt5Agg as a driver for matplotlib, please install it.

		PyQt5: https://pypi.org/project/PyQt5/

3. Copy and modify example data and configuration files in the working
   directory

   Data file: A comma- or tab-separated data sheet has at least six
      columns of 207Pb/235U, 207Pb/235U_error, 206Pb/238U,
      206Pb/238U_error, 207Pb/206Pb, and 207Pb/206Pb_error

   Configuration file: A text file describing variables used by this
      script.  Th file name is assumed to be the same with that of the
      data file except the extension which is cfg.

[Usage]:

After installation of libraries listed above, you can run the script
in a terminal window, for examples,

	$ python3 UPbplot.py -i data.csv -c data.cfg
	$ python3 UPbplot.py -i data.csv -d pdf -f
	$ python3 UPbplot.py -d qt5agg

The script assumes the configuration file name is "data.cfg" as
defaults, if the input data file name is "data.csv".

Command-line options:

	Options:
	-h, --help                  Show this help message and exit
	-i FILE, --in=FILE          Name of input data file
	-c FILE, --cfg=FILE         Name of configure file
	-o FILE, --out=FILE         Name of output file (when pdf driver is used)
	-d DRIVER, --driver=DRIVER  Choose from [pdf (default), Qt5Agg, TKAgg, macosx]
	-f, --force-overwrite       Force overwrite the pre-existing pdf

[Features]:

1. There are some manners to separate discordant data from concordant ones.  This script can automatically classify data whose error ellipses intersect the concordia curves as discordant data.  Set ``disc_type: 5'' in the configuration file to use this method.

2. Generalized ESD test is a common procedure to identify outliers in a data set.  This script can automatically exclude outliers at a certain critical value of the confidince level (eg., "outlier_alpha: 0.05").  Set "opt_outlier: True" to use this option.

3. Correction for initial disequilibria can be applied based on Sakata (2017, 2018), when ``opt_correction_disequilibrium: True''

4. Correction for common Pb by 207Pb method can be applicable, when ``opt_correction_common_Pb: True''

5. Upper and lower intercept ages can be calculated.


[References]:

Sakata et al., 2017, A new approach for constraining the magnitude of initial disequilibrium in {Quaternary} zircons by coupled uranium and thorium decay series dating. Quaternary Geochronology, vol. 38, p. 1--12, https://doi.org/10.1016/j.quageo.2016.11.002.

Sakata, S., 2018, A practical method for calculating the U-Pb age of Quaternary zircon: Correction for common Pb and initial disequilibria. Geochemical Journal, vol. 52, p. 281--286, https://doi.org/10.2343/geochemj.2.0508.


