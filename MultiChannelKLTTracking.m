function [u,v,label,metric] = MultiChannelKLTTracking(multiChannelImages,detectedPoints,weight,maxIterations,convergence,windowSize,numPyramid,maxDisplacement,outputMetric)
%MULTICHANNELKLTTRACKING Summary of this function goes here
%   Detailed explanation goes here

% mulitChannelImages : nRows * nCols * nChannels * nFrames

% check input

img = double(multiChannelImages); % save a copy in the scope

imageDimensions = size(img);
nRows = imageDimensions(1); nCols = imageDimensions(2);
if length(imageDimensions) > 3
    nChannels =  imageDimensions(3);
    nFrames = imageDimensions(4);
elseif length(imageDimensions) == 3
    nChannels = 1;
    nFrames = imageDimensions(3);
else
    error("input images are neither multichannel, nor 2-frame")
end

if nFrames ~= 2
    error("input images have more than 2 frames, if tracking sequential frames" + ...
        "is requried, please use this function iteratively over every 2-frame stack");
end



if size(detectedPoints,2) == 2
    label = ones(size(detectedPoints,1));
elseif size(detectedPoints,2) == 3
    label = detectedPoints(:,3);
else
    error(strcat("the second input arguement, detectedPoints, has wrong dimensions.",...
        "it should be either n-by-2, i.e., [x,y], or n-by-3, i.e., [x,y,l] ", ...
        "where x and y are coordinates, l is the label, having the same length of x and y, that indicates the original channel of a detected point"));
end

nPoints = size(detectedPoints,1);


% by default the first frame is the reference frame

x_pos = detectedPoints(:,1);
y_pos = detectedPoints(:,2);

guess_lvl = zeros(length(x_pos),2);
for iPyramid = numPyramid:-1:0
    % resize images for each pyramid level
    scale = 0.5^iPyramid;
    img_lvl = imresize(img,scale,'bilinear','Antialiasing',true);
    x_lvl = x_pos * scale; y_lvl = y_pos * scale;
    [u_lvl,v_lvl] = multiChannelKLTTracking_base();
    if iPyramid ~= 0
    guess_lvl = 2*[u_lvl, v_lvl];
    end
end

u = u_lvl;
v = v_lvl;

if nargin == 8
    outputMetric = true;
end

if outputMetric == true
metric = normxcorr2metric();
end


    function [u,v] = multiChannelKLTTracking_base()
        imref = img_lvl(:,:,:,1);
        imdef = img_lvl(:,:,:,2);
        nRows_lvl = size(imref,1);
        nCols_lvl = size(imref,2);
        halfSize = round(windowSize/2);
        x_lvl_new = zeros(nPoints,1);
        y_lvl_new = zeros(nPoints,1);
        parfor i = 1:nPoints
            x = x_lvl(i);
            y = y_lvl(i);

            [ymesh,xmesh] = ndgrid(max(1, y - halfSize):min(nRows_lvl, y + halfSize),max(1, x - halfSize):min(nCols_lvl, x + halfSize));
            y1 = max(1,ceil(min(min(ymesh)))); y2 = min(ceil(max(max(ymesh))),nRows_lvl);
            x1 = max(1,ceil(min(min(xmesh)))); x2 = min(ceil(max(max(xmesh))),nCols_lvl);

            Ix = [];
            Iy = [];
            I0 = [];

            for ich = 1:nChannels
                imref_window_current_ch = imref(y1:y2,x1:x2,ich);
                [Gx,Gy] = getGradient(imref_window_current_ch);
                [ymesh_grad,xmesh_grad] = ndgrid(y1:y2,x1:x2);
                gx = griddedInterpolant(ymesh_grad,xmesh_grad,Gx,'cubic','cubic');
                gy = griddedInterpolant(ymesh_grad,xmesh_grad,Gy,'cubic','cubic');
                Iold = griddedInterpolant(ymesh_grad,xmesh_grad,imref_window_current_ch,'cubic','cubic');
                Ix = [Ix;weight(i,ich)*gx(ymesh,xmesh)];
                Iy = [Iy;weight(i,ich)*gy(ymesh,xmesh)];
                I0 = [I0;weight(i,ich)*Iold(ymesh,xmesh)];
            end

            A = [sum(Ix.*Ix,"all") sum(Ix.*Iy,"all"); sum(Ix.*Iy,"all") sum(Iy.*Iy,"all")];

            delta = guess_lvl(i,:);
            for iter = 1:maxIterations+1
                % update displaced mesh
                xmesh = xmesh + delta(1);
                ymesh = ymesh + delta(2);
                x = x + delta(1);
                y = y + delta(2);
                y1 = max(1,ceil(min(min(ymesh)))); y2 = min(ceil(max(max(ymesh))),nRows_lvl);
                x1 = max(1,ceil(min(min(xmesh)))); x2 = min(ceil(max(max(xmesh))),nCols_lvl);
                [ymesh_grad,xmesh_grad] = ndgrid(y1:y2,x1:x2);
                if (norm(delta) > 0 && norm(delta) < convergence) || isempty(xmesh_grad) || ...
                        isempty(ymesh_grad) || max(y1,y2) < 1 ||...
                        min(y1,y2) > nRows ||max(x1,x2) < 1 ...
                        || min(x1,x2) > nCols || ...
                        x1 == x2||y1==y2 || iter > maxIterations || norm([x-x_lvl(i) y-y_lvl(i)]) > maxDisplacement
                    break;
                end

                I1 = [];
                for ich=1:nChannels
                    imdef_window_current_ch = imdef(y1:y2,x1:x2,ich);
                    Inew = griddedInterpolant(ymesh_grad,xmesh_grad,imdef_window_current_ch,'cubic','cubic');
                    I1 = [I1; weight(i,ich) * Inew(ymesh,xmesh)];
                end
                It = I0-I1;
                b = [sum(It.*Ix,"all"); sum(It.*Iy,"all")];
                delta = pinv(A)*b;
            end
            x_lvl_new(i) = x;
            y_lvl_new(i) = y;
        end
        u = x_lvl_new - x_lvl;
        v = y_lvl_new - y_lvl;
    end

    function metric = normxcorr2metric()
        imref = img(:,:,:,1);
        imdef = img(:,:,:,2);
        halfSize = round(windowSize/2);
        metric = struct("spatialDeviation",[],"maximumCorrelation",[]);
        for ich = 1:nChannels
            metric(ich) = struct("spatialDeviation",zeros(nPoints,1),"maximumCorrelation",zeros(nPoints,1));
            imref_current_ch = imref(:,:,ich);
            imdef_current_ch = imdef(:,:,ich);
            parfor i = 1:nPoints
                xrange_ref = [max(1, ceil(x_pos(i) - halfSize)), min(nCols, ceil(x_pos(i) + halfSize))];
                yrange_ref = [max(1, ceil(y_pos(i) - halfSize)), min(nRows, ceil(y_pos(i) + halfSize))];
                xrange_def = [max(1, ceil(x_pos(i) + u(i) - halfSize)), min(nCols, ceil(x_pos(i) + u(i) + halfSize))];
                yrange_def = [max(1, ceil(y_pos(i) + v(i) - halfSize)), min(nRows, ceil(y_pos(i) + v(i) + halfSize))];
                templateImPatch = imref_current_ch(yrange_ref(1):yrange_ref(2),xrange_ref(1):xrange_ref(2));
                offsetImPatch = imdef_current_ch(yrange_def(1):yrange_def(2),xrange_def(1):xrange_def(2));
                
                [xpeak,ypeak,corrpeak] = xcorrPeakGaussFit(templateImPatch,offsetImPatch,7);
                spatialDeviation(i) = (xpeak^2 + ypeak^2)^0.5;
                maximumCorrelation(i) = corrpeak;
            end
            metric(ich).spatialDeviation = spatialDeviation';
            metric(ich).maximumCorrelation = maximumCorrelation';
        end
    end

end




