clc
clear all
close all

% *************************************************************************
% Lucas & Kanade (1981)
%   Estimate and display optic flow for an expansion motion.
%
%   Copyright (C) 2013  Florian Raudies, 01/14/2013, Boston University.
%   License, GNU GPL, free software, without any warranty.
% *************************************************************************

% Generate stimulus.
addpath('./../'); % add root directory to path to access auxiliary functions.
opt.h       = 128; % set parameters for stimulus.
opt.w       = 128;
opt.fNum    = 15;
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
len1    = sample/maxSpeed1;
len2    = sample/maxSpeed2;

% *************************************************************************
% Estimate and display flow for the 1st sequence.
% *************************************************************************
opt.DiffT = reshape([-1 0 1],[1 1 3]);
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
figure('Position',[150 50 600 600]);
quiver(X(IndexY,IndexX),        Y(IndexY,IndexX),...
       Dx2(IndexY,IndexX)*len2, Dy2(IndexY,IndexX)*len2,0,'-k');
axis equal ij; axis([-10 w+10 -10 h+10]);
title(sprintf(['Expansion.\n',...
    'Esimated optic flow, sampled %d times for a %d x %d',...
    'pixels resolution.'], sample,h,w));
drawnow;
