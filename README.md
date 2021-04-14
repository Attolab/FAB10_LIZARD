# FAB10_LIZARD
Implementation of the LIZARD technique on FAB10 beamline

## Installation

This was intended to be a pymodaq pid model. But as we started the development with the version 1.6.3 of pymodaq ([http://pymodaq.cnrs.fr/en/latest/index.html](url)) we fixed some bugs quite in a nasty way and commented some lines in the pymodaq repository, in particular in the pid_controller.py file (and made also modifications in the pymodaq_plugins repository). The good way would have been to to modify only the pymodaq_pid_models repository. This will be done in a later version.
So in this repository we include the pymodaq_pid_models repository but also pymodaq repository and the pymodaq plugins as they were when we tested the all feedback system.

Probably the best way to install from scratch (we mainly follow the pymodaq installation procedure [http://pymodaq.cnrs.fr/en/latest/usage/Installation.html](url)):

Install Miniconda or Anaconda and cd to the *condabin* folder.

Create and activate a new python 3.7 environment

    conda create -n new_environment_name python=3.7
    conda activate new_environment_name

Install pyqt with conda

    conda install pyqt

Install pymodaq version 1.6.3 with pip

    pip install pymodaq==1.6.3

Copy paste in the *site-packages* folder of your environment the folders that are contained in this repository.


## System requirements

Windows 10 (for now the used plugins are not compatible with Linux, but pymodaq should be)
python 3.7
pymodaq 1.6.3
pymodaq plugins : daq_move_SmarActMCS and daq_1Dviewer_LecroyWaverunner6Zi

Delay stage: Smaract SLC with S option (sensor module with nanometer precision is absolutely necessary) and MCS controller (MCS2 will probably not be compatible as it is)
Oscilloscope: Lecroy Waverunner 6Zi

## References

[http://pymodaq.cnrs.fr/en/latest/index.html](url)
[https://journals.aps.org/prapplied/abstract/10.1103/PhysRevApplied.15.034036](url)
[https://aip.scitation.org/doi/10.1063/5.0032116](url)
