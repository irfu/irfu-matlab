irfu-matlab [![Build Status](https://travis-ci.com/irfu/irfu-matlab.svg?branch=master)](https://travis-ci.com/irfu/irfu-matlab)
===========

Acknowledgment
--------------

Please use the following to acknowledge use of IRFU-Matlab in your publications:

> Data analysis was performed using the IRFU-Matlab analysis package available at https://github.com/irfu/irfu-matlab 

Examples
--------

A few examples of irfu-matlab usage can be seen under web page https://sites.google.com/site/irfumatlab/


Installation
-------------

Download the code and add to the Matlab path the directory 'irfu-matlab'.

To download the code from the command line run

```sh
git clone git://github.com/irfu/irfu-matlab.git
```
If you want only the latest version and not the full repository run
```sh
git clone --depth=1 git://github.com/irfu/irfu-matlab.git
```
Any time you want to update the code to the latest version, run from the command line 
```sh
git pull
```
If you do not use command line there are github programs for Windows and Mac, see github web page. 

Usage
-----

Each time starting new Matlab session execute in Matlab:

```matlab
irf
```

To see help execute in Matlab 
```matlab
help irfu-matlab
```

To test an example execute in Matlab:
```matlab
Example_1
```
this should download Cluster data from the CSA archive and plot the figure consisting of 5 panels with data from different instruments. 

Issues
-----

If you experience any issues with irfu-matlab please [submit an issue](https://github.com/irfu/irfu-matlab/issues) via our github page.

When submitting a new issue, please provide information required to replicate the issue as well as information regarding operating system and version of Matlab used.
