clc
clear all
close all

% *************************************************************************
% Horn & Schunck (1981)
%   Estimate the optic flow for the Yosemite with clouds image sequence.
%   The partial differential equation (PDE) system is solved by using the
%   successive over relaxation (SOR) method, an iterative method. The 
%   displayed estimated optic flow is updated during these iterations.
%
%   Copyright (C) 2013  Florian Raudies, 01/14/2013, Boston University.
%   License, GNU GPL, free software, without any warranty.
% *************************************************************************

% Load Yosemite with clouds image sequence.
addpath('./../');
filePattern = './../YosemiteWithClouds/ImgFrame%05d.pgm';
ImgSeq = readImgSeq(filePattern,0,1);
rmpath('./../');


% Estimate optic flow for the image sequence.
opt.eta = 0.1;
[Dx Dy] = estimateOpticFlow2D(ImgSeq,opt);

% Display the estimated optic flow.
h       = size(ImgSeq,1);
w       = size(ImgSeq,2);
[Y X]   = ndgrid(1:h, 1:w); % pixel coordinates.
sample  = 8;
IndexX  = 1:sample:w;
IndexY  = 1:sample:h; 
len     = sample*2;
figure('Position',[50 50 600 600]); 
quiver(X(IndexY,IndexX),      Y(IndexY,IndexX),...
       Dx(IndexY,IndexX)*len, Dy(IndexY,IndexX)*len,0,'-k');
axis equal ij; axis([-10 w+10 -10 h+10]);
title(sprintf(['Esimated optic flow, sampled %d times for a %d x %d',...
    ' pixels resolution.'], sample,h,w));
