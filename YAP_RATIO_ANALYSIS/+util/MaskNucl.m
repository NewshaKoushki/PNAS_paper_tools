%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name of file: MaskNucl
% Version 1.0

% Purpose: to segment a Nucleius in the image using thresholding 

% Input variables: 
%   im: the image
% 
%   output variable: 
%   imout: the mask

% created: Feb 20, 2017 12:07 PM
% Author: Rosa Kaviani
% Revisions: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [imout]= MaskNucl (imin)

imin = medfilt2(imin,[7 7]);
[~, threshold] = edge(imin, 'Prewitt');
fudgeFactor =1;
BWs = edge(imin,'sobel', threshold * fudgeFactor);
se90 = strel('line', 3, 90);
se0 = strel('line', 3, 0);
BWsdil = imdilate(BWs, [se90 se0]);


BWdfill = imfill(BWsdil, 'holes');

BWnobord = imclearborder(BWdfill, 4);

seD = strel('diamond',1);
BWfinal = imerode(BWnobord,seD);
BWfinal = imerode(BWfinal,seD);

imout = bwareaopen(BWfinal,250);
end