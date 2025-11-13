%%
% 
% * This example script demonstrates how to use our Multi-channel KLT tracking * 
%  
% 

% clear workspace and display history
clear
clc
% load example data from a simulated traction field
load("example_data.mat","multiChannelImages");
% mulitChannelImages should be a 4D array( nRows * nCols * nChannels *
% nFrames) containing pixels and frames
%% DETECTION
% define a detection method
% detectionCase 1：minimum eigen feature detection 
%               2: local maximum finder
% note that both detection methods require MATLAB Computer Vision Toolbox

detectionCase = 1;

switch detectionCase
    case 1 
    % parameteres work for the example data.
        detectionMethod = "MinEigenFeatures";
        MinQuality = 0.01;
        FilterSize = 5;
        [x,y,labels] = PointDetection(multiChannelImages,detectionMethod,MinQuality,FilterSize);
    case 2
    % Did not test these parameters for example data
        detectionMethod = "LocalMaixmaFinder";
        maxnum = 1e5;
        neighborsize = [5,5];
        threshold = 30.0;
        [x,y,labels] = PointDetection(multiChannelImages,detectionMethod,maxnum,neighborsize,threshold);
end


% output structure of PointDetection:
% x,y are the spatial coordiantes of detected points.
% label indicates the channel a point is detected from.


%% TRACKING

% define the parameters for multi-channel Lucas-Kanade method
nPts = length(x); % get the number of points 
nChs = size(multiChannelImages,3); % get the number of channels

detectedPoints = [x,y,labels]; % pack input positions
weight = ones(nPts,nChs); % initialize the weighting matrix, ones(nPts,nChs) for unweighted tracking

maxIterations = 30; % the internal iteration ends until the maximum number of iterations is reached
convergence = 0.01; % the internal iteraction ends until the convergence is reached (pix) 
windowSize = 32; % interogation window size (pix)
numPyramid = 2; % number of pyramid level, 0 -> original size, n-> (1/2)^n of original size
maxDisplacement = 2; % maximum displacement magnitude (pix) at each pyramid level

% select if iterative weighting scheme is used, weighting is given by the
% maximum of normalized cross-correlation of the interogation window of
% each detected point.
iterativeWeighting = 1;
%%
switch iterativeWeighting
    case 0
    [u,v,labels,metric] = MultiChannelKLTTracking(multiChannelImages,detectedPoints,weight,maxIterations,convergence,windowSize,numPyramid,maxDisplacement);
    case 1
        nWeightedIterations = 2; % the number of weighted iterations
        [u,v,labels,metric] = MultiChannelKLTTrackingWeightedIteration(multiChannelImages,detectedPoints,nWeightedIterations,maxIterations,convergence,windowSize,numPyramid,maxDisplacement);
end

% tracking output strcutures:
% u,v are the displacement components in x and y direction, respectively.
% You can use quiver(x,y,u,v) to visualize the displacement vector plot.
% labels indicates the original channel for the points.
% metric characterizes the deviation beteween the maximum cross-correlation
% peaks and the tracked positions.
% metric has a length of the number of channels. For each channel, it has
% a field, spatialDeviation and maxCorrelation. You can use customized
% thresholds to exculde poorly tracked vectors.




