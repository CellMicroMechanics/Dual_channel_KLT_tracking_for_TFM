function [x,y,labels] = PointDetection(multiChannelImages,Method,varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% Currently only works for 2-frame images
% mulitChannelImages : nRows * nCols * nChannels * nFrames
imageDimensions = size(multiChannelImages);
if length(imageDimensions) > 2
    nChannels = imageDimensions(3);
else
    nChannels = 1;
end
% nRows = imageDimensions(1);
% nCols = imageDimensions(2);

img = multiChannelImages;
x = [];
y = [];
labels = [];
switch Method
    case "LocalMaximaFinder"
        switch length(varargin)
            case 0
                maxnum = 1e5;
                neighborsize = [5,5];
                threshold = 30;
            case 1
                maxnum = varargin{1};
                neighborsize = [5,5];
                threshold = 30;
            case 2
                maxnum = varargin{1};
                neighborsize = varargin{2};
                threshold = 30;
            case 3
                maxnum = varargin{1};
                neighborsize = varargin{2};
                threshold = varargin{3};
        end     
        for i = 1:nChannels
            LocalMaximumFinder_tmp = vision.LocalMaximaFinder("MaximumNumLocalMaxima",maxnum,"NeighborhoodSize",neighborsize,"Threshold",threshold);
            ptstmp = LocalMaximumFinder_tmp(img(:,:,i,1));
            release(LocalMaximumFinder_tmp);
            x = [x;ptstmp(:,1)];
            y = [y;ptstmp(:,2)];
            labels = [labels;i*ones(size(ptstmp.Location),1)];
            %clear ptstmp;
        end
    case "MinEigenFeatures"
        switch length(varargin)
            case 0
                MinQuality = 0.01;
                FilterSize = 5;
            case 1
                MinQuality = varargin{1};
                FilterSize = 5;
            case 2
                MinQuality = varargin{1};
                FilterSize = varargin{2};
        end
        for i = 1:nChannels
            ptstmp = detectMinEigenFeatures(img(:,:,i,1),"MinQuality",MinQuality,"FilterSize",FilterSize);
            x = [x;ptstmp.Location(:,1)];
            y = [y;ptstmp.Location(:,2)];
            labels = [labels;i*ones(size(ptstmp.Location,1),1)];
            %clear ptstmp;
        end
end

x = double(x);
y = double(y);


end