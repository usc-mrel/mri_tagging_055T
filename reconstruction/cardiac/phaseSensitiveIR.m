function [ imgCorrected, imgPhase, imgWindowed ] = phaseSensitiveIR( images, winSize )
% phaseSensitiveIR: Corrects the phase of the image using the phase of a
% low resolution image by taking a central region of k-space.
% The central region is defined by the winSize. The phase of the low resolution image
% is then used to correct the phase of the entire image.
%
%   INPUTS:
%       images:     3D array of images
%       winSize:    Size of the window used to correct the phase
%
%   OUTPUTS:
%       imgCorrected:   Corrected image
%       imgPhase:       Phase of the low resolution image
%       imgWindowed:    Windowed image

if nargin<2
    winSize = 4;
end
[imgSize,~,~] = size(images);

xMin = floor(imgSize/2-winSize/2+1);
yMin = floor(imgSize/2-winSize/2+1);

kspace =  fftshift(fft2(images));
kspaceCtr = kspace(xMin:xMin+winSize-1, yMin:yMin+winSize-1, :);

imageLR = ifft2(fftshift(padarray(kspaceCtr,[(imgSize-winSize)/2, (imgSize-winSize)/2],0,'both')));

imgPhase = angle(imageLR);
imgCorrected = images./exp(1i*imgPhase);
imgWindowed = imgCorrected;

imgCorrected = real(imgCorrected);
end
