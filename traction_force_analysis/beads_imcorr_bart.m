function [x, y, dx, dy] = beads_imcorr_bart(im1, im2)

% BEADS_IMCORR calculates the folder_for_displacement field between two images.
%	The program is optimized for images of fluorescent beads in the gel.
%	Input:	images im1 and im2, blocksize_start (size of the largest
%				correlation window), resolution (size of the smallest
%				correlation window).
%	Output:	x (x-coordinates), y (y-coordinates),
%				dx (x-displacements), dy (y-displacements),
%				saved in a file in the folder e:\Cells\Displacements.
%				The file name contains the name of the current folder
%				and names of the two images.

%	Iva Marija Tolic-Norrelykke 03-21-01

%======================================================================================
% PARAMETERS (feel free to change them!)

blocksize_start	    =  32;		% Resolution at the starting level; 32
resolution       	=  16;		% Resolution at the final level; Changed from 8 to 16
min_pixval			=	50;		% Smallest pixel value for a bead
min_corr			=	0.5;		% Threshold value of the cross-correlation
method				=	'v4';		% Method for interpolation
maxNsteps			=	3;			% Max # of steps in searching for the matching window

%======================================================================================
% INITIAL p1, dx, blocksize

p1 = im1;
p2 = im2;
%if whichcell==4			%Apr-06-01\cell4, crop is not a square, add zeros to it
%    m1 = size(p1,1);
%    n1 = size(p1,2);
%   if m1>n1
%      p1(:,(n1+1):m1)	=	0;
%  elseif m1<n1
%      p1((m1+1):n1,:)=	0;
%  end
%
%   	m2 = size(p2,1);
%    n2 = size(p2,2);
%   if m2>n2
%      p2(:,(n2+1):m2)	=	0;
%  elseif m2<n2
%     p2((m2+1):n2,:)=	0;
%  end
% %end    						% (whichcell)
% im1=p1;
% im2=p2;

dx = zeros(round(size(p1)./blocksize_start));
dy = dx;

blocksize = blocksize_start;

%======================================================================================
% PROCESS p1, p2 ON BLOCKS

while blocksize >= resolution,
    
    [x,y] = meshgrid((blocksize/2 + 0.5) : blocksize : (size(im1,2)-blocksize/2 + 0.5),...
        (blocksize/2 + 0.5) : blocksize : (size(im1,1)-blocksize/2 + 0.5));
    
    step_number = 1;
    suspect_i 	= NaN;
    clear bx by;
    
    while ~isempty(suspect_i) && step_number <= maxNsteps,
        
        % bxn:  NEW DISPLACEMENTS FROM BLOCKS
        [bxn,byn] = disp_on_blocks(p1,p2,blocksize,0);
        
        % bx: FINAL DISPLACEMENTS AT THIS LEVEL (blocksize), NOT JUST THE RESIDUAL ONES
        if ~exist('bx', 'var'),
            bx = zeros(size(bxn)) + round(dx);
            by = zeros(size(byn)) + round(dy);
        end;
        
        % bxa: NEW DISPLACEMENTS TO BE ADDED TO bx
        bxa = zeros(size(bx));
        bya = zeros(size(by));
        
        if isfinite(suspect_i),
            for i = 1 : length(suspect_i),
                bxa(suspect_i(i), suspect_j(i)) = bxn(i);
                bya(suspect_i(i), suspect_j(i)) = byn(i);
                bxa = reshape(bxa,size(bx));
                bya = reshape(bya,size(bx));
            end; %(for i)
        else
            bxa = bxn;
            bya = byn;
        end;
        
        % ADDING OLD AND NEW DISPLACEMENTS
        bx = bx + bxa;
        by = by + bya;
        
        % FINDING NEW DISPLACEMENTS ~= 0 AND MOVING ONLY THOSE BLOCKS
        [suspect_i, suspect_j] = find(round(bxa) ~= 0 | round(bya) ~= 0);
        [p1, p2, bx, by] = move_blocks(im1,im2,bx,by,suspect_i,suspect_j);
        
        % ENDING THE LOOP: IF step_number >= maxNsteps OR DISPLACEMENTS CONVERGED
        if step_number > 2,
            if (bx == bxoldold) & (by == byoldold), step_number = 1000; end;
        end; %(if step_number)
        if step_number > 1,
            if (bx == bxold) & (by == byold), step_number = 100; end;
            bxoldold = bxold;
            byoldold = byold;
        end; %(if step_number)
        bxold = bx;
        byold = by;
        step_number = step_number + 1;
        
    end; %(while ~isempty(suspect_i) & step_number < maxNsteps)
    
    %======================================================================================
    % CHECK POINT
    
    [p1, p2, bx, by] = move_blocks(im1,im2,bx,by);
    p1 = blkproc(p1,[size(im1,1) size(p1,2)], @(x) x(:));
    p1 = reshape(p1,size(im1));
    p2 = blkproc(p2,[size(im2,1) size(p2,2)], @(x) x(:));
    p2 = reshape(p2,size(im2));
    
    [bxn,byn,cor_whole] = disp_on_blocks(p1,p2,blocksize,1);
    cor_whole = blkproc(cor_whole,[size(im1,1) size(cor_whole,2)], @(x) x(:));
    cor_whole = reshape(cor_whole,size(im1));
    max_wh 	 = blkproc(cor_whole,[blocksize blocksize],'max(x(:))');
    
    max_pix1 = blkproc(p1,[blocksize blocksize],'max(x(:))');
    
    bad = find(round(bxn) ~= 0 | round(byn) ~= 0 | ...  				% not best alligned
        abs(bx-dx)	> 	blocksize/2 | abs(by-dy) > blocksize/2 | ...	% moved too much
        max_pix1		< 	min_pixval | ...	  									% too black
        max_wh		< 	min_corr); 												% max corr too low
    
    bx = bx + bxn;
    by = by + byn;
    
    bx(bad) = NaN;
    by(bad) = NaN;
    
    %======================================================================================
    
    if blocksize > resolution,
        blocksize = blocksize/2;
    end;
    
    %======================================================================================
    % MAKE dx AND dy FOR THE NEXT LEVEL
    
    [xn,yn] = meshgrid((blocksize/2 + 0.5) : blocksize : (size(im1,2)-blocksize/2 + 0.5),...
        (blocksize/2 + 0.5) : blocksize : (size(im1,1)-blocksize/2 + 0.5));
    
    fin = find(isfinite(bx));
    
    dx = griddata(x(fin),y(fin),bx(fin),xn,yn,method);
    dy = griddata(x(fin),y(fin),by(fin),xn,yn,method);
    
    
    dx(~isfinite(dx)) = 0;
    dy(~isfinite(dy)) = 0;
    
    %======================================================================================
    % MAKE p1, p2 FOR THE NEXT LEVEL
    
    if size(x,1) ~= size(xn,1),
        [p1, p2, dx, dy] = move_blocks(im1,im2,round(dx),round(dy));
        p1 = blkproc(p1,[size(im1,1) size(p1,2)], @(x) x(:));
        p1 = reshape(p1,size(im1));
        p2 = blkproc(p2,[size(im2,1) size(p2,2)], @(x) x(:));
        p2 = reshape(p2,size(im2));
    else
        blocksize = 0;
    end;
    
end; %(while blocksize >= resolution)

%Commented on 7/22/08 in order to compute the tensile strain fields.
%Removed the mean correction
%dx = dx - mean2(dx);
%dy = dy - mean2(dy);

dx = dx;
dy = dy;

%======================================================================================
% TAKE THE EVEN NUMBER OF POINTS ON THE GRID

s1 = size(x,1);
s2 = size(x,2);
if mod(s1,2) ~= 0,
    x  = x (1:s1-1,:);
    y  = y (1:s1-1,:);
    dx = dx(1:s1-1,:);
    dy = dy(1:s1-1,:);
end

if mod(s2,2) ~= 0,
    x  = x (:,1:s2-1);
    y  = y (:,1:s2-1);
    dx = dx(:,1:s2-1);
    dy = dy(:,1:s2-1);
end

end


function [p1new,p2new,dx_after,dy_after] = move_blocks(p1,p2,dx,dy,suspect_i,suspect_j)

% MOVE_BLOCKS creates images p1new and p2new from the images p1 and p2
%	such that at each position in p2new there is a block from p2 that is the
%	the best match to the block in p1new at the same position.

%	Iva Marija Tolic-Norrelykke 03-21-01

%======================================================================================

if nargin <= 4,
    suspect_i = [];
    suspect_j = [];
end;

dx(find(~isfinite(dx))) = 0;
dy(find(~isfinite(dy))) = 0;
dx = round(dx);
dy = round(dy);

blocksize = size(p2,1) / size(dy,1);

[xwhole,ywhole] = meshgrid(1 : size(p2,2), 1 : size(p2,1));

%======================================================================================
% MOVE BLOCKS IN p2

x_beg = blkproc(xwhole,[blocksize blocksize],'min(x(:))');
x_end = blkproc(xwhole,[blocksize blocksize],'max(x(:))');
y_beg = blkproc(ywhole,[blocksize blocksize],'min(x(:))');
y_end = blkproc(ywhole,[blocksize blocksize],'max(x(:))');

% CORNERS IN p1
corners_old(:,:,1) = y_beg;
corners_old(:,:,2) = y_end;
corners_old(:,:,3) = x_beg;
corners_old(:,:,4) = x_end;

% CORNERS IN p2 TO COMPARE WITH p1
corners(:,:,1) = y_beg + dy;
corners(:,:,2) = y_end + dy;
corners(:,:,3) = x_beg + dx;
corners(:,:,4) = x_end + dx;

% CORNERS THAT GO OUTSIDE OF p2 HAVE TO GET INSIDE
[fc1i, fc1j] = find(corners(:,:,1) < 1);
corners(fc1i, fc1j, 1) = 1;
corners(fc1i, fc1j, 2) = blocksize;
[fc3i, fc3j] = find(corners(:,:,3) < 1);
corners(fc3i, fc3j, 3) = 1;
corners(fc3i, fc3j, 4) = blocksize;

[fc2i, fc2j] = find(corners(:,:,2) > size(p2,1));
corners(fc2i, fc2j, 1) = corners(fc2i, fc2j, 1) - (corners(fc2i, fc2j, 2) - size(p2,1));
corners(fc2i, fc2j, 2) = size(p2,1);
[fc4i, fc4j] = find(corners(:,:,4) > size(p2,2));
corners(fc4i, fc4j, 3) = corners(fc4i, fc4j, 3) - (corners(fc4i, fc4j, 4) - size(p2,2));
corners(fc4i, fc4j, 4) = size(p2,2);

% FINAL CORNERS IN p2
[in, jn, kn] = find(~isfinite(corners));
corners(in,jn,:) = corners_old(in,jn,:);
corners_new = corners;

%======================================================================================
% CALCULATE dx_after, dy_after

dx_after = corners_new(:,:,3) - corners_old(:,:,3);
dy_after = corners_new(:,:,1) - corners_old(:,:,1);

%======================================================================================
% MAKE NEW p2

newpict = [];
oldpict = [];

if ~isempty(suspect_i),
    corners_old1 = [];
    corners_new1 = [];
    for i = 1 : length(suspect_i),
        corners_old1 = [corners_old1; corners_old(suspect_i(i),suspect_j(i),:)];
        corners_new1 = [corners_new1; corners_new(suspect_i(i),suspect_j(i),:)];
    end; %(for i)
    corners_old = corners_old1;
    corners_new = corners_new1;
end; %(isempty)

for j = 1 : size(corners_new,2),
    for i = 1 : size(corners_new,1),
        oldpict = [oldpict; p1(corners_old(i,j,1):corners_old(i,j,2),...
            corners_old(i,j,3):corners_old(i,j,4))];
        newpict = [newpict; p2(corners_new(i,j,1):corners_new(i,j,2),...
            corners_new(i,j,3):corners_new(i,j,4))];
    end; %(for i)
end; %(for j)

p1new = oldpict;
p2new = newpict;

end