# LaueG
Archive of code used to process data from the Neutron Laue Instrument, KOALA, at the Australian Centre for Neutron Scattering

The Australian Nuclear Science and Technology Organisation decided to defund the KOALA instrument along with the two scientists that developed and operated the instrument. Rumours are that funds may be obtained sometime in the future to restart KOALA, though the original two scientists will no longer be involved. I created this repository to archive the software suite, developed almost entirely by myself, which processes the raw data from the KOALA instrument into standard *.hkl data files as used by crystallography refinement packages SHELX, CRYSTALS, JANA, etc.

The code is highly specific and is only applicable to Laue single-crystal neutron-diffraction data from imaging-plate or other "storage" detectors. In particular, it is unsuitable for time-of-flight data from a spallation source. The code runs under MS Windows, though emulation on other systems are said to work.

The GUI and graphical display of the data, along with some minor calculations, are performed using the open-source program, SCILAB. Complex calculations are performed by single instances of CMD running individual Fortran executables. Calculation options and input/output data are transferred between SCILAB and the Fortran executables using temporary files of text and binary data. SCILAB is free software and contains several major faults that have not been addressed since reported (even when fixes are given) for more than 10 years. I implemented these fixes myself and compiled several Java modules to overwrite the standard ones during the installation of LaueG. For this reason LaueG only works for one specific version of SCILAB, v5.5.2. The Fortran code was compiled using Intel Fortran Compiler XE 13.1 under MS Visual Studio 2010. The code is written in F77 using Tab-Line format for the source code. I emphasised robustness in the algorithms with the hope that 90% of the data sets can be processed using default options, 9% require tweaking of parameters as possible by a skilled operator, and 1% need my own expertise and possible code changes. I seem to have got those percentages about right.

```
Literature:
"Accurate data processing for neutron Laue diffractometers"
  R. O. Piltz, J. Appl. Cryst. (2018). 51, 635-645.
"LaueG software for displaying and processing neutron Laue images"
  R. O. Piltz, J. Appl. Cryst. (2018). 51, 963-965.
The Laue spot integration software is a highly-modified version of ARGONNE_BOXES:
  Wilkinson, C., Khamis, H. W., Stansfield, R. F. D. & McIntyre, G. J. (1988). J. Appl. Cryst. 21, 471–478.
NLSCON is used as the nonlinear least-squares solver:
  "Newton Methods for Nonlinear Problems: Affine Invariance and Adaptive Algorithms" ,
  Deuflhard, P. (2004), Heidelberg: Springer-Verlag.
```
A user manual for the LaueG suite is included in the documentation section of this archive.

==================

I am a scientist not a software engineer. I had to "reinvent the wheel" in terms of managing this package as I did not know about such things at the time. This software was written out of necessity and has been under continuous development since its inception. I know that certain functions have not been completed and do not work properly. I did my best!

Finally, the name LaueG was inspired by Ali G, a fictional character created by the comedian Sacha Baron Cohen.

Dr Ross Piltz, March 2026
