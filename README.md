irfu-matlab
===========

Examples
--------

A few examples of irfu-matlab usage can be seen under web page https://sites.google.com/site/irfumatlab/


Installation
-------------

Download the code and add to the Matlab path the directory 'irfu-matlab'.

To download the code from the command line run

> git clone git://github.com/irfu/irfu-matlab.git

If you want only the latest version and not the full repository run

> git clone --depth=1 git://github.com/irfu/irfu-matlab.git

Any time you want to update the code to the latest version, run from the command line 

> git pull

If you do not use command line there are github programs for Windows and Mac, see github web page. 

Usage
-----

Each time starting new Matlab session execute in Matlab:

> irf

To see help execute in Matlab 

> help irfu-matlab

To test an example execute in Matlab:

> Example_1

this should download Cluster data from the CSA archive and plot the figure consisting of 5 panels with data from different instruments. 

Issues
-----

If you experience any issues with irfu-matlab please [submit an issue](https://github.com/irfu/irfu-matlab/issues) via our github page.

When submitting a new issue, please provide information required to replicate the issue as well as information regarding operating system and version of Matlab used.
