# UPbplot.py

[Introduction]

This is a script for calculation and visualization tool of U-Pb age
data.  It enables us to

1. plot scattered data with error ellipses on conventional
(207Pb/235U-206Pb/238U) and Terra--Wasserburg (238U/206Pb-207Pb/206Pb)
concordia diagrams, and

2. calculate one- or two-dimensional weighted mean, concordia, and
concordia-intercept ages with errors on both concordia diagrams.


[License, Author, Copyright, Publication date]

License: Apache License, Version 2.0
Author: Atsushi Noda
Copyright: Geological Survey of Japan, AIST
Publication date: Jul 03, 2019 (ver 0.0.9)


[Citation]

Noda, Atsushi (2016) UPbplot.py: A python script for U-Pb age data
analysis. Open-File Report, no. 634, Geological Survey of Japan, AIST.

Noda, Atsushi (2017) A new tool for calculation and visualization
of U-Pb age data: UPbplot.py.  Bulletin of the Geological Survey of
Japan, 68(3), 131-140.


[Development environment]

The script was written in Python 3.6.8 and tested on MacOSX
(Mojave).  If you use it under Python 2, please modify lines using 
``print'' function in the script.  Comments are added to the relevant lines.

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

$ python UPbplot.py -n -i data.csv -c data.cfg
$ python UPbplot.py -n -i data.csv -d pdf -f
$ python UPbplot.py -d qt5agg


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

