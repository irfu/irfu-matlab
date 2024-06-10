# irfu-matlab [![ci-build](https://github.com/irfu/irfu-matlab/actions/workflows/ci-build.yml/badge.svg)](https://github.com/irfu/irfu-matlab/actions/workflows/ci-build.yml)  [![DOI](https://zenodo.org/badge/8636532.svg)](https://zenodo.org/doi/10.5281/zenodo.11550090)

## Acknowledgment

Please use the following to acknowledge use of IRFU-Matlab in your publications:

> Data analysis was performed using the IRFU-Matlab analysis package available at https://github.com/irfu/irfu-matlab

## Examples

A few examples of irfu-matlab usage can be seen under web page https://sites.google.com/site/irfumatlab/

## Installation

Download the code and add to the Matlab path the directory 'irfu-matlab'. Two methods are described below, a plain basic install and a more development oriented install (for contributors).

#### Basic install

If you intend to use irfu-matlab without pushing code then to simply download the code from the command line run

```sh
git clone https://github.com/irfu/irfu-matlab.git
```

If you want only the latest version and not the full repository run

```sh
git clone --depth=1 https://github.com/irfu/irfu-matlab.git
```

Any time you want to update the code to the latest version, run from the command line

```sh
git pull
```

If you do not use command line there are github programs for Windows and Mac, see github web page.

#### Development install

This section is mostly for those of you who intend to develop and contribute code yourself (i.e. `git push`) and not just use 'irfu-matlab'.

It is recommended you first add SSH keys to your github.com account, please see github web pages for details on how to add an SSH key to your account.
And then, when SSH keys have been configured, instead of cloning using the basic `https://` command above instead from the command line run

```sh
git clone ssh://git@github.com/irfu/irfu-matlab.git
```

Then finally before actually contributing code, please have a look at our [contributing](https://github.com/irfu/irfu-matlab/blob/master/.github/CONTRIBUTING.md) guidelines.

## Usage

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

## Issues

If you experience any issues with irfu-matlab please [submit an issue](https://github.com/irfu/irfu-matlab/issues) via our github page.

When submitting a new issue, please provide information required to replicate the issue as well as information regarding operating system and version of Matlab used.
