function d = track_optical_flow(img1,img2,minQuality,MaxErr,BlkSz,pos)

if nargin == 5
    pts = detectMinEigenFeatures(img1,"MinQuality",minQuality);
    pos = pts.Location;
elseif nargin == 6
    clear minQuality
end


tracker = vision.PointTracker("MaxBidirectionalError",MaxErr,"BlockSize",BlkSz);
initialize(tracker,pos,img1)
[pos2,valid] = tracker(img2);
release(tracker)

d = double([pos(valid,:),pos2(valid,:)-pos(valid,:)]);

figure;
im(:,:,1) = img1; im(:,:,2) = img2;
sliceViewer(im);
hold on;
quiver(d(:,1),d(:,2),d(:,3),d(:,4),0,'r');

end