# Code for: Dataset of human-single neuron activity during a Sternberg working memory task 

[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
[![Generic badge](https://img.shields.io/badge/release-1.0.0-green.svg)](https://github.com/rutishauserlab/workingmem-release-NWB/releases/1.0.0)
[![Generic badge](https://img.shields.io/badge/DOI-Status_Pending-orange.svg)](https://www.youtube.com/watch?v=MOu4_wpy40U)

## Introduction

This repository contains tools/methods to import, format, and analyze a sample dataset of human single-neuron electrophysiological and behavioral data. This dataset is formatted entirely in [Neurodata Without Borders (NWB)](https://www.nwb.org/) format. This format has extensive [documentation](https://nwb-schema.readthedocs.io/en/latest/index.html) and support in both MATLAB and Python environments. 

This code accompanies the paper: (Citation Pending). [Link to paper](https://www.youtube.com/watch?v=MOu4_wpy40U)

Abstract of the paper:
>We present a dataset of 1809 single neurons recorded from the human medial temporal lobe (amygdala and hippocampus) and medial frontal lobe (anterior cingulate cortex, pre-supplementary motor area, ventral medial prefrontal cortex) across 41 sessions from 21 patients that underwent intracranial monitoring for epileptic activity. Subjects first performed a screening task (907 neurons), based on which we identified images for which highly selective cells were present in the medial temporal lobe. Subjects then performed a working memory task (902 neurons), in which they were sequentially presented with 1-3 images, and following a maintenance period, were asked if a probe was identical to one of the currently maintained images. This Neurodata Without Borders (NWB) formatted dataset includes spike times, extracellular spike waveforms, stimuli presented, behavior, electrode locations, and subject demographics. As validation, we replicate previous findings on the existence of concept cells and their persistent activity during working memory maintenance.  This dataset provides a substantial amount of rare human single neuron recordings together with behavior, thereby enabling investigation of the neural mechanisms of working memory at the single-neuron level.


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

`dandi version 0.55.1, nwbinspector version 0.4.28, & pynwb 2.3.1.` 

[PyNWB Documentation](https://pynwb.readthedocs.io/en/stable/)



## MATLAB Analysis








## Contributors
* [Michael Kyzar](mailto:kyzarnexus@gmail.com)
* [Ueli Rutishauser](mailto:Ueli.Rutishauser@cshs.org) (Principal Investigator)

>([Citation Pending](https://www.youtube.com/watch?v=MOu4_wpy40U))

## Funding

Acquisition of this dataset has been supported by the National Institute of Mental Health (grants U01NS117839, U01NS098961, and R01MH110831)

## License 

"workingmem-release-NWB" Copyright (c) 2019, Rutishauser Lab. All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


