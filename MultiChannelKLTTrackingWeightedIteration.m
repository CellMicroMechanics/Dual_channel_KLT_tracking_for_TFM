function [u,v,label,metric] = MultiChannelKLTTrackingWeightedIteration(multiChannelImages,detectedPoints,nWeightedIterations,maxIterations,convergence,windowSize,numPyramid,maxDisplacement)

if numel(size(multiChannelImages))> 3
nChs = size(multiChannelImages,3);
elseif numel(size(multiChannelImages)) == 3
    nChs = 1;
end
outputMetric = true;
weight = ones(size(detectedPoints,1),nChs);
if nChs > 1

    for i = 1:nWeightedIterations
        [u,v,label,metric] = MultiChannelKLTTracking(multiChannelImages,detectedPoints,weight,maxIterations,convergence,windowSize,numPyramid,maxDisplacement,outputMetric);
        weight = horzcat(metric.maximumCorrelation);
    end
else
     [u,v,label,metric] = MultiChannelKLTTracking(multiChannelImages,detectedPoints,weight,maxIterations,convergence,windowSize,numPyramid,maxDisplacement,outputMetric);
end


end

