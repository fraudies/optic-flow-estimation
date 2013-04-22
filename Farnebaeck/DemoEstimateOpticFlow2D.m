clc
clear all
close all

% *************************************************************************
% Farnebaeck (2000).
% This is a re-implementation of code provided by Gunnar Farnebaeck, which 
% had been avaiable at http://www.isy.liu.se/~gf/software/ some years ago.
% Unlike to Farnebaeck's implmentation, I do not use mex files.
% 
%   Copyright (C) 2013  Florian Raudies, 01/14/2013, Boston University.
%   License, GNU GPL, free software, without any warranty.
% *************************************************************************

% Load the Yosemite with coulds sequence.
addpath('./../'); % add root directory to path to access auxiliary functions.
ImgSeq = readImgSeq('./../YosemiteWithClouds/ImgFrame%05d.pgm',0,14);
maxSpeed = 5; % Set to a reasonable value, might not be correct.
rmpath('./../');

% The default parameters assume that intensities range from 0 to 255.
[Dx Dy] = estimateOpticFlow2D(ImgSeq*255);

% Display the estimated optic flow.
h       = size(ImgSeq,1);
w       = size(ImgSeq,2);
[Y X]   = ndgrid(1:h, 1:w); % pixel coordinates.
sample  = 8;
IndexX  = 1:sample:w;
IndexY  = 1:sample:h;
% For the display the flow is scaled by division with its maximum speed and
% multiplication with the sampling factor.
len     = sample/maxSpeed;

figure('Position',[50 50 600 600]);
quiver(X(IndexY,IndexX),      Y(IndexY,IndexX),...
       Dx(IndexY,IndexX)*len, Dy(IndexY,IndexX)*len,0,'-k');
axis equal ij; axis([-10 w+10 -10 h+10]);
title(sprintf(['Yosemite sequence with clouds.\n',...
    'Esimated optic flow, sampled %d times for a %d x %d',...
    'pixels resolution.'], sample,h,w));
