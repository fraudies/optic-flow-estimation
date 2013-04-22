clc
clear all
close all

% *************************************************************************
% Adelson & Bergen (1985).
%   Estimation of optic flow for rotation and expansion motion.
%
%   Copyright (C) 2013  Florian Raudies, 01/08/2013, Boston University.
%   License, GNU GPL, free software, without any warranty.
% *************************************************************************

% Generate the first two image sequences.
addpath('./../'); % add root directory to path to access auxiliary functions.
opt.h       = 128; % Set parameters for stimulus.
opt.w       = 128;
opt.fNum    = 19;
opt.omega   = pi/300;
[ImgSeq1 maxSpeed1] = randomTextureRotation(opt);
opt.expRate = 1.005;
[ImgSeq2 maxSpeed2] = randomTextureExpansion(opt);
rmpath('./../');

% Display the estimated optic flow.
h       = size(ImgSeq1,1);
w       = size(ImgSeq1,2);
[Y X]   = ndgrid(1:h, 1:w); % pixel coordinates.
sample  = 8;
IndexX  = 1:sample:w;
IndexY  = 1:sample:h;
% For the display the flow is scaled by division with its maximum speed and
% multiplication with the sampling factor.
len1    = sample/maxSpeed1;
len2    = sample/maxSpeed2;


% *************************************************************************
% Estimate and display flow for the 1st sequence.
% *************************************************************************
[Dx1 Dy1] = estimateOpticFlow2D(ImgSeq1);
figure('Position',[50 50 600 600]);
quiver(X(IndexY,IndexX),        Y(IndexY,IndexX),...
       Dx1(IndexY,IndexX)*len1, Dy1(IndexY,IndexX)*len1,0,'-k');
axis equal ij; axis([-10 w+10 -10 h+10]);
title(sprintf(['Rotation.\n',...
    'Esimated optic flow, sampled %d times for a %d x %d',...
    'pixels resolution.'], sample,h,w));
drawnow;

% *************************************************************************
% Estimate and display flow for the 2nd sequence.
% *************************************************************************
[Dx2 Dy2] = estimateOpticFlow2D(ImgSeq2);
figure('Position',[350 50 600 600]);
quiver(X(IndexY,IndexX),        Y(IndexY,IndexX),...
       Dx2(IndexY,IndexX)*len2, Dy2(IndexY,IndexX)*len2,0,'-k');
axis equal ij; axis([-10 w+10 -10 h+10]);
title(sprintf(['Expansion.\n',...
    'Esimated optic flow, sampled %d times for a %d x %d',...
    'pixels resolution.'], sample,h,w));
