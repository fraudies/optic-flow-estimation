clc
clear all
close all

% *************************************************************************
% Otte & Nagel (1994).
%   Estimate the optic flow for the Yosemite with clouds sequence.
%   Copyright (C) 2013  Florian Raudies, 01/14/2013, Boston University.
%   License, GNU GPL, free software, without any warranty.
% *************************************************************************

% Load Yosemite with clouds.
addpath('./../');
ImgSeq = readImgSeq('./../YosemiteWithClouds/ImgFrame%05d.pgm',0,14);
maxSpeed = 5; % Set to a reasonable value, which might not be the true max.
rmpath('./../');

% Estimate optic flow for the image sequence.
[Dx Dy] = estimateOpticFlow2D(ImgSeq);

% Display the estimated optic flow.
h       = size(ImgSeq,1);
w       = size(ImgSeq,2);
[Y X]   = ndgrid(1:h, 1:w); % pixel coordinates.
sample  = 8;
IndexX  = 1:sample:w;
IndexY  = 1:sample:h; 
len     = sample/maxSpeed;

figure('Position',[50 50 600 600]); 
quiver(X(IndexY,IndexX),      Y(IndexY,IndexX),...
       Dx(IndexY,IndexX)*len, Dy(IndexY,IndexX)*len,0,'-k');
axis equal ij; axis([-10 w+10 -10 h+10]);
title(sprintf(['Esimated optic flow, sampled %d times for a %d x %d',...
    'pixels resolution.'], sample,h,w));

