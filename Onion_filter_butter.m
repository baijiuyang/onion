function [vOutput, outFiltered, outExtended] = Carrot4_filter_butter (Hz,vInput,cutoff,order,pad)
% filter_butter   Butterworth low-pass filter with linear extrapolation. 
%   filter_butter(hZ, varargin) Filters and returns input vector with a 
%   Butterworth low-pass filter using standard VENLab parameters and 
%   sample rate hZ. Extrapolates first and last 0.5 sec to avoid transient 
%   at either end. 

% Filter both x and y independently.

% Created by Henry Harrison, 2010
% Modified by Kevin Rio, 2013
% Modified by Gregory Dachner, 2017 / gdachner@gmail.com
% Modified by Jiuyang Bai, 2017
% Brown University

% Compute parameters

% standard VENLab parameters is cutoff = 0.6


[B,A] = butter(order,cutoff/(Hz/2));

% Pad data with timesteps to extrapolate into
tot_time = 1:length(vInput)+pad*2;

extrap_pre = 1:pad;
extrap_post = length(tot_time)-pad+1:length(tot_time);

exist_data_location = length(extrap_pre)+1:extrap_post(1)-1;
exist_data_location = exist_data_location'; % ' means transpose

toExtrap = vertcat(extrap_pre',extrap_post');

% Extrapolate
outExtrap = interp1(exist_data_location,vInput,toExtrap,'linear','extrap');
outExtended = vertcat(outExtrap(1:pad),vInput,outExtrap(pad+1:end));

% Apply butterworth filter and return input to original size
outFiltered = filtfilt(B,A,outExtended); 
vOutput = outFiltered(exist_data_location);

