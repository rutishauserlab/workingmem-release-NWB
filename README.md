# Source code for: Dataset of human-single neuron activity during a Sternberg working memory task 

[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
[![Generic badge](https://img.shields.io/badge/release-1.0.0-green.svg)](https://github.com/rutishauserlab/workingmem-release-NWB/releases/1.0.0)
[![Generic badge](https://img.shields.io/badge/DOI-Status_Pending-orange.svg)](https://youtu.be/MOu4_wpy40U)

## Introduction

This repository contains the code that accompanies our data release Kyzar et al. 'Dataset of human-single neuron activity during a Sternberg working memory task'. The purpose of the code in this repository is to provide examples of how to use the released data. This dataset is formatted in the [Neurodata Without Borders (NWB)](https://www.nwb.org/) format, which can easily be accessed from both MATLAB and Python as described [here](https://nwb-schema.readthedocs.io/en/latest/index.html) . 

This code accompanies the following data descriptor: (Citation Pending). [Link to paper(Pending)](https://youtu.be/MOu4_wpy40U) 

The papers that describe the scientific results that are based on this dataset are:
(add here Kaminski 2017, 2020)

Abstract of the paper:
>We present a dataset of 1809 single neurons recorded from the human medial temporal lobe (amygdala and hippocampus) and medial frontal lobe (anterior cingulate cortex, pre-supplementary motor area, ventral medial prefrontal cortex) across 41 sessions from 21 patients that underwent intracranial monitoring for epileptic activity. Subjects first performed a screening task (907 neurons), based on which we identified images for which highly selective cells were present in the medial temporal lobe. Subjects then performed a working memory task (902 neurons), in which they were sequentially presented with 1-3 images, and following a maintenance period, were asked if a probe was identical to one of the currently maintained images. This Neurodata Without Borders (NWB) formatted dataset includes spike times, extracellular spike waveforms, stimuli presented, behavior, electrode locations, and subject demographics. As validation, we replicate previous findings on the existence of concept cells and their persistent activity during working memory maintenance.  This dataset provides a substantial amount of rare human single neuron recordings together with behavior, thereby enabling investigation of the neural mechanisms of working memory at the single-neuron level.

<p align="center">
  <img width="400" height="500" src="https://github.com/rutishauserlab/workingmem-release-NWB/blob/main/assets/git_task_overview_white.png">
  <img width="400" height="500" src="https://github.com/rutishauserlab/workingmem-release-NWB/blob/main/assets/git_brain_areas.png">
</p>


## Installation (Code)

This repository can be downloaded by entering the following commands:

`cd $target_directory`

`git clone https://github.com/rutishauserlab/workingmem-release-NWB.git`

## Installation (MatNWB)

Running the provided code and analyzing the dataset in MATLAB requires the download and initialization of MatNWB, a MATLAB interface for reading and writing NWB 2.x files. Instructions for how to [download and initialize MatNWB](https://github.com/NeurodataWithoutBorders/matnwb) have been listed on the project's public git repo. Further documentation for how to use MatNWB can be found [here](https://neurodatawithoutborders.github.io/matnwb/). MatNWB version 2.6.0.2 was used for the curation and analysis of this dataset.

## Installation (Data)

The dataset is available in NWB format from the Dandi Archive, under [Dandiset #469](https://dandiarchive.org/dandiset/000469).


Dandi datasets are accessible through the Dandi command line interface (CLI). To install this Python client, use `pip install dandi` or `conda install -c conda-forge dandi`, depending on your Python environment setup. 


After installing the Dandi CLI, use `dandi download https://dandiarchive.org/dandiset/000469` to download the dataset. 

## Installation & File Validation (Python)

NWB Files can additionally be loaded and analyzed using the [PyNWB](https://github.com/NeurodataWithoutBorders/pynwb) python package. Further documentation can be found [here](https://pynwb.readthedocs.io/en/stable/). 


Validation of this dataset was performed using PyNWB (2.3.1) and PyNWB-dependent packages, such as nwbinspector (0.4.28) and dandi (0.55.1). The command lines used for each method are as follows:
* dandi: `dandi validate $target_directory`
* nwbinspector: `nwbinspector $target_directory`
* PyNWB: `Get-ChildItem $target_directory -Filter *.nwb -Recurse | % { $_.FullName }; python -m pynwb.validate $file_list`
  <!--- Test the PyNWB method again. There seems to be an access error --->

All validators returned no errors in data formatting & best-use practices across all uploaded files. 


## MATLAB Analysis

The main script in this repo, `NWB_SB_import_main.m`, is designed to analyze the released dataset and to reproduce the figures & metrics noted in Kyzar et. al. 2023, Kaminski et. al. 2017, and Kaminski et. al. 2020. It can calculate several metrics related to behavior (reaction time, accuracy), spike sorting, and single-unit (SU) activity during the screening (SC) & Sternberg (SB) tasks.

### Steps to Use the Script
* **Set Parameters:** The first section of the script sets important parameters. The `taskFlag` is used to specify which tasks are accessed and are defined as `1` (SCREENING), `2` (STERNBERG), or `3` (BOTH). The `importRange` is the range of subject IDs for the dataset. For the current release, subject IDs have a range of `1:21`. 

* **Initialization and Pathing:** The script then defines the directory paths for the code, the currently installed MatNWB package, and the dataset, and then adds them to the MATLAB path. If figures are generated, there is an additional option to add a custom save destination. Please ensure that the defined paths in the script are correct for your setup. This section also uses MatNWB's generateCore() function to initialize the NWB API if it has not been initialized already.

* **Import Datasets From Folder:** The script will then import datasets from the given folder using the `NWB_importFromFolder` function. Only files specified using `importRange` and `taskFlag` will be loaded into the workspace. 

* **Extracting Single Units:** Single unit information is extracted from the loaded NWB files for ease of indexing, using the `NWB_SB_extractUnits` function. If spike waveforms are not needed for analysis, the `load_all_waveforms` flag can be set to `0` to only extract the mean waveform. All future plots will use this mean waveform instead of a spike waveform pdf. 

* **Screening Analysis:** This section is preceded by a parameters section, which allows for the control of various stages of the analysis and plotting process. For example, one can choose to plot figures for significant cells by setting `paramsSC.doPlot = 0` or filter units being used for analysis by specifying a minimum firing rate threshold `paramsSC.rateFilter`. To disable screening analysis of all cells entirely, set `paramsSC.calcSelective = 0`. 

* **Sternberg Analysis:**: This section is also preceded by a parameters section, with additional controls included for the type of cell type to plot (`1`: Concept, `2`: Maint, `3`: Probe, `4`: All). Similar to screening, analysis can be disabled by setting `paramsSB.calcSelective = 0`.

* **Example Neurons:** Additional sections have been added for each task type that optionally plots the example cells found in Kyzar et al 2023 & Kaminski et al 2017. To reduce load times for NWB files, set `importRange = [4 7 15 16 21]` for screening and `importRange = [7 14 16]` for Sternberg.

* **Selectivy by Area:** The script also calculates the proportion of selective cells by area and plots bar/pie charts comparing the screening and Sternberg tasks. This can be disabled by setting `plotAreas = 0`.

* **Spike Sorting Quality Metrics:** This section plots spike sorting metrics for single units recorded in the Sternberg/screening tasks. These metrics include the percentage of inter-spike intervals (ISIs) that were less than 3 ms, mean firing rates for all units, coefficient of variation (CV2) values for all units, signal-to-noise ratio (SNR) of the peak of the mean waveform, mean SNR for all values in a unit’s mean waveform, pairwise projection distance between each unit in which multiple units were found on the same electrode, isolation distance (scaled to log 10 for ease of viewing) across all units for which this metric was defined.


Please make sure to thoroughly read the comments in the code to understand the functionality of each part. If you encounter any problems, please report them as issues in the repository.


## Contributors
* [Michael Kyzar](mailto:kyzarnexus@gmail.com)
* [Ueli Rutishauser](mailto:Ueli.Rutishauser@cshs.org) (Principal Investigator)

>([Citation Pending](https://youtu.be/MOu4_wpy40U))

## Funding

Acquisition of this dataset has been supported by the National Institute of Mental Health (grants U01NS117839, U01NS098961, and R01MH110831)

## License 

"workingmem-release-NWB" Copyright (c) 2019, Rutishauser Lab. All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


