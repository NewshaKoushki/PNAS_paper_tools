function [bx,by,corsh] = disp_on_blocks(im1,im2,blockSize,subpix)

% DISP_ON_BLOCKS calculates x- and y-displacements for each distinct block
%	taken from the image im1. The size of the blocks is determind by the
%	parameter blocksize. For each block in the image im1 the corresponding
%	block in the image im2 is found (i.e., the block at the same position),
%	and the cross-correlation function between the two blocks is formed.
%	Coordinates of the peak of the cross-correlation function constitute
%	the displacement vector of the block. (Displacements go from im1 to im2.)
%	The parameter subpix determines whether the resulting values of displacements
%	are integers or not; if subpix > 0, the values of displacements are
%	non-integers, otherwise they are integers.

%	bx is a matrix of x-displacements, by a matrix of y-displacements, and
%	corsh the fft-shifted cross-correlation by block between im1 and im2.

%	Iva Marija Tolic-Norrelykke 03-21-01

%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

blocksize = [blockSize blockSize];

if nargin < 4,
    subpix	= 	0;
end;

% SUBTRACT THE MEAN OF EACH BLOCK
im1 = blkproc(im1, blocksize, @(x) x - mean2(x));
im2 = blkproc(im2, blocksize, @(x) x - mean2(x));

%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% CALCULATE THE CROSS-CORRELATION FUNCTION

autocor_whole1 = blkproc(im1, blocksize, @(x) ifft2((fft2(x)).*conj(fft2(x))));
autocor_whole2 = blkproc(im2, blocksize, @(x) ifft2((fft2(x)).*conj(fft2(x))));
autocor_max1   = blkproc(autocor_whole1, blocksize, @(x)  repmat(x(1,1),size(x)));
autocor_max2   = blkproc(autocor_whole2, blocksize, @(x) repmat(x(1,1),size(x)));
cor_whole1     = blkproc(im1, blocksize, @(x) conj(fft2(x)));
cor_whole2     = blkproc(im2, blocksize, @(x) fft2(x));
cor_whole      = cor_whole2 .* cor_whole1;
cor_whole		= blkproc(cor_whole, blocksize, @(x) ifft2(x));
cor_whole      = real(cor_whole) ./ sqrt(autocor_max1.*autocor_max2);
cor_whole(~isfinite(cor_whole)) = zeros(size(find(~isfinite(cor_whole))));
corsh          = blkproc(cor_whole, blocksize, @(x) fftshift(x));

%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% FIND THE PEAK OF THE CROSS-CORRELATION FUNCTION

if subpix > 0,
    bx				= blkproc(corsh, blocksize, @(x) center_x_1d(x));
    by				= blkproc(corsh, blocksize, @(x) center_y_1d(x));
else
    bx				= blkproc(corsh, blocksize, @(x) max_cor_bx(x));
    by				= blkproc(corsh, blocksize, @(x) max_cor_by(x));
end; %(if subpix > 0)

% blkproc, blockproc

end

%%
function  bx = max_cor_bx(x)

[~, bx] = find(x == max(x(:)));
bx = bx(1) - size(x,2)/2 - 1;

if length(bx) > 1
   bx = 0;
end

end

%%
function  by = max_cor_by(x)

[by, ~] = find(x == max(x(:)));
by = by(1) - size(x,1)/2 - 1;

if length(by) > 1
   by = 0;
end

end


%%
function [cx] = center_x_1d(x)

% CENTER_X_1D calculates the x-position of the peak of the 
%	cross-correlation function x by fitting the nearest 
% 	neighbors (left and right) of the maximum value of x to a parabola.

try
   
[mi,mj]	= find(x == max(max(x)));
mi			= mi(1);
mj			= mj(1); 

[xobj, yobj]	= 	meshgrid(1 : size(x,2), 1 : size(x,1));

x_fit = xobj(mi, mj-1 : mj+1);
y_fit = x(mi, mj-1 : mj+1);

[a,b] = polyfit(x_fit,y_fit,2);

x_vals = [min(x_fit) : 0.01 : max(x_fit)];
y_vals = polyval(a, x_vals);

cx = x_vals(find(y_vals == max(y_vals))) - size(x,2)/2 - 1;

if length(cx) > 1, cx = mean(cx); end

catch
   cx = NaN;
end
end

%%
function [cy] = center_y_1d(x)

% CENTER_Y_1D calculates the y-position of the peak of the 
%	cross-correlation function x by fitting the nearest 
% 	neighbors (up and down) of the maximum value of x to a parabola.

try
   
[mi,mj]	= find(x == max(max(x)));
mi			= mi(1);
mj			= mj(1);

[xobj, yobj]	= 	meshgrid(1 : size(x,2), 1 : size(x,1));

x_fit = yobj(mi-1 : mi+1, mj);
y_fit = x(mi-1 : mi+1, mj);

[a,b] = polyfit(x_fit,y_fit,2);

x_vals = [min(x_fit) : 0.01 : max(x_fit)];
y_vals = polyval(a, x_vals);

cy = x_vals(find(y_vals == max(y_vals))) - size(x,1)/2 - 1;

if length(cy) > 1, cy = mean(cy); end

catch
   cy = NaN;
end

end