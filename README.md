# Millimetro
 Sample data and processing backend for Millimetro, a long-range tag-based mmWave localization system.

Millimetro is a tag-based localization system which operates accurately at long ranges to accurately localize mulitple modulating tags. This repository contains the MATLAB files required to run it as well as sample collected data files to demonstrate its capabilities. The sample_code.m file contains the signal processing for the back-end of Millimetro with other .m files used as helper functions. The data folder has 4 different datasets: a static dataset with long-range tags, a dataset with a mobile radar and static tags, a dataset with tags moving across angle of arrival, and a dataset with two tags while colocated in range operating at different modulation frequencies.

This code has been released online for transparency and in order to encourage further development in the area.

## Installation/Requirements

The Millimetro code is run on MATLAB, version R2021a or later (earlier installations may work). To install, download all files included with this repository into a folder and open sample_code.m in MATLAB. Set the downloaded folder as your working directory, then run the code. 

## Features
Running the code will take all frames of data in some collection and perform a set of operations in this order: a set of Fourier Transforms, matched filtering to locate the tag's range,  matched filtering to locate the tag's angle of approach, and then comparisons to ground truth data, if it was collected. It is recommend to put a breakpoint at the end of the code in order to visualize the results from each data collection.

Listed below are some of the parameters of interest.
| Variable | Index | Meaning |
| ------ | ------ | ----- |
| experiments | 1:4 | Changes which data files the code indexes into |
| jj | - | Changes which file within a set of experiments the code processes |
| rangeAll | (jj, frame number, tag number (if multiple)) | Calculated range of the tag, if located |
| angleAll | (jj, frame number, tag number (if multiple)) | Calculated angle of the tag, if located |
| tagRangeErr | (jj, frame number, tag number (if multiple)) | Calculated error of the tag to gt |
| tagAoAErr | (jj, frame number, tag number (if multiple)) | Calculated error of the tag to gt |
| RPExtAll | (range bin, chirp number bin, frame number, antenna) | The range-profile data from the radar |
| RDAll | (range bin, doppler bin, frame number, antenna) | The range-doppler data from the radar |

The code is currently configured to process all of the data from a file at once in order to improve performance. It is fairly easy to restructure it from its current state into a version which runs frame by frame while collecting data (i.e. the data frames are processed using a causal system).

##File List
| File name | Description |
| ------ | ----- |
| sample_code | Main file to handle millimetro processing |
| importCSV | Imports data from a .csv file |
| findPeaks2D | Helper function to locate the peak of the matched filter |
| aruco_calc |  Derives grund truth parameters from imported .csv file |
| tagRangeErr | Calculated error of the tag to gt |
| tagAoAErr | Calculated error of the tag to gt |
| RPExtAll | The range-profile data from the radar |
| RDAll | The range-doppler data from the radar |
