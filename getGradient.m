function [Gx,Gy] = getGradient(image,sigma,method) 

if nargin == 1
sigma = 1.0;
method = "sobel";
elseif nargin == 2
    method = "sobel";
end

[Gx,Gy] = imgradientxy(imgaussfilt(image,sigma),method);

end

