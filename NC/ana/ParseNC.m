%% ParseNC.m
% Summary
%   Goes through NC data, extracts and calculates relevant information, and
% saves the results. Can run each section separately after running 'SET
% VARIABLES' and 'LOAD META DATA.' Note some sections require data from a
% previous section. This is an example script of how to use the functions.
%
% 1. Sleep
%   Extracts time stamp of sleep using delta oscillations (saved as .mat
% file) and accelerometer signals. Saved as Sleep.mat.
%
% 2. Sorting
%   Finds parameters for sorting for each recording channel and displays the
% obtained parameters. Manual input required for each channel on keeping or
% discarding the spike. Assumes there is only one spike per channel. After
% obtaining the parameters, it sorts for each channel detected and saves
% the spike times in separate .mat files.
%
% 3. stLFP
%   Plots spike triggered LFPs over the entire duration of the experiment for
% each spike and for each recording channel.
%
% RJY 06/22/2018

%% Set variables
clear; close all;
path = 'Y:\~NeuroWest\Spanky\Neurochip';
day = 'S_20180613_07';

%% Load meta data
[fpath,fname,Channels,fs,session_time] = getNCData(path,day);

%% Do functions
% Sleep
saveSleep(fpath,fname,session_time,28);

% Sorting 
saveSortingParams(fpath,fname,Channels,fs,600);
saveSortedSpikes(fpath,fname,session_time,fs);

% Plot spike stats
plotSpikeStats(fpath);

% stLFP
saveSTLFP(fpath,fname,day,session_time,Channels,fs,[15,50]);





