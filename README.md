# Millimetro
 Sample data and processing backend for Millimetro, a long-range tag-based mmWave localization system

Millimetro is a tag-based localization system which operates accurately at ranges up to 50 m. This repository contains the MATLAB files required to run it as well as a collection of sample collected data files to demonstrate it's capabilities.

## Installation/Requirements

The Millimetro code is run on MATLAB, version R2021a or later (earlier installations may work). To install, download all files in this repository, and open up pipeline_sample.m in MATLAB. Set the downloaded folder as your working directory, the run the code. 

## Features
Running the code will take all frames of data and perform a set of operations in this order: a set of FFTs, matched filtering to locate the tag's range,  matched filtering to locate the tag's angle of approach, and then comparisons to ground truth data, if it was collected. We recommend putting a breakpoint at the end of the code in order to visualize the results from each collection of frames.
Listed below are some of the parameters of interest.
| Variable | Index | Meaning |
| ------ | ------ | ----- |
| rangeAll | (jj, frame number, tag number (if multiple)) | |
| angleAll | (jj, frame number, tag number (if multiple)) | |
| jj | - | Changes which file within a set of experiments the code processes |
| experiments | 1:4 | Changes which data files the code indexes into |
| RPExtAll | (range bin, chirp number bin, frame number, antenna) | The range-profile data from the radar |
| RDAll | (range bin, doppler bin, frame number, antenna) | The range-doppler data from the radar |