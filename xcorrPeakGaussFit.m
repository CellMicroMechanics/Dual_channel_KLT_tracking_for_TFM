function [xpeak,ypeak,varargout] = xcorrPeakGaussFit(img1,img2,fitRadius)
    
    if ~isequal(size(img1),size(img2))&& ~ismember(1,size(img1)) &&  ~ismember(1,size(img2)) && ~isempty(img1) && ~isempty(img2)
        [nr1,nc1] = size(img1);
        [nr2,nc2] = size(img2);
        [rs1,cs1] = ndgrid(1:nr1,1:nc1);
        [rs2,cs2] = ndgrid(1:nr2,1:nc2);
        [rs,cs] = ndgrid(1:min(nr1,nr2),1:min(nc1,nc2));
        tmp1 = griddedInterpolant(rs1,cs1,img1,'linear','linear');
        tmp2 = griddedInterpolant(rs2,cs2,img2,'linear','linear');
        img1tmp = tmp1(rs,cs);
        img2tmp = tmp2(rs,cs);
    elseif ismember(1,size(img1)) ||  ismember(1,size(img2)) || isempty(img1) || isempty(img2)
        xpeak = 0;
        ypeak = 0;
        for i = 1:nargout-2
            varargout{i} = 0;
        end
        return;
    else
        img1tmp = img1;
        img2tmp = img2;
    end

    
 if length(unique(img1tmp)) > 1 
    cc1 = normxcorr2(img1tmp,img2tmp);
    [r1max,c1max] = find(cc1 == max(cc1(:)),1);
    rlow1 = max(1,r1max - fitRadius);
    rhigh1 = min(size(cc1,1),r1max + fitRadius);
    fitRangeY = rlow1:rhigh1;
    
    clow1 = max(1,c1max - fitRadius);
    chigh1 = min(size(cc1,2),c1max + fitRadius);

    fitRangeX = clow1:chigh1;

    cc = cc1(fitRangeY,fitRangeX);
    
    [xtmp,ytmp] = meshgrid(fitRangeX,fitRangeY); 

    peakFitted = lsqcurvefit(@Gauss2D,[1,1,c1max,r1max],[xtmp(:) ytmp(:)],cc(:),[],[],optimset("Display","off"));

    xpeak = peakFitted(3) - size(img1,2);
    ypeak = peakFitted(4) - size(img1,1);
 else
     ypeak = 0;
     xpeak = 0;
     peakFitted = zeros(1,5);
 end
    
    switch nargout - 2
        case 1
            varargout{1} = peakFitted(1);
        case 2
            varargout{1} = peakFitted(1); % AMPLITUDE
            varargout{2} = peakFitted(2); % SIGMA
        case 3
            varargout{1} = peakFitted(1);
            varargout{2} = peakFitted(2);
            varargout{3} = peakFitted(5); % residual
    end
end