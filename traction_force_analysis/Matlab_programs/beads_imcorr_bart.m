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

dx = zeros(size(p1,1)/blocksize_start, size(p1,2)/blocksize_start);
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
   
	while ~isempty(suspect_i) & step_number <= maxNsteps,
                    
      % bxn:  NEW DISPLACEMENTS FROM BLOCKS
      [bxn,byn] = disp_on_blocks(p1,p2,blocksize,0);
      
      % bx: FINAL DISPLACEMENTS AT THIS LEVEL (blocksize), NOT JUST THE RESIDUAL ONES
      if ~exist('bx'),
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
         if bx == bxoldold & by == byoldold, step_number = 1000; end;
      end; %(if step_number)
      if step_number > 1,
         if bx == bxold & by == byold, step_number = 100; end;
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
   p1 = blkproc(p1,[size(im1,1) size(p1,2)],'mat2vertcol(x)');
   p1 = reshape(p1,size(im1));
   p2 = blkproc(p2,[size(im2,1) size(p2,2)],'mat2vertcol(x)');
   p2 = reshape(p2,size(im2));

   [bxn,byn,cor_whole] = disp_on_blocks(p1,p2,blocksize,1);
   cor_whole = blkproc(cor_whole,[size(im1,1) size(cor_whole,2)],'mat2vertcol(x)');
   cor_whole = reshape(cor_whole,size(im1));
   max_wh 	 = blkproc(cor_whole,[blocksize blocksize],'max2(x)');
   
   max_pix1 = blkproc(p1,[blocksize blocksize],'max2(x)');
   
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


	dx(find(~isfinite(dx))) = 0;
	dy(find(~isfinite(dy))) = 0;
   
%======================================================================================
% MAKE p1, p2 FOR THE NEXT LEVEL

	if size(x,1) ~= size(xn,1),
		[p1, p2, dx, dy] = move_blocks(im1,im2,round(dx),round(dy));
   	p1 = blkproc(p1,[size(im1,1) size(p1,2)],'mat2vertcol(x)');
   	p1 = reshape(p1,size(im1));
   	p2 = blkproc(p2,[size(im2,1) size(p2,2)],'mat2vertcol(x)');
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
end;
if mod(s2,2) ~= 0,
   x  = x (:,1:s2-1);
   y  = y (:,1:s2-1);
   dx = dx(:,1:s2-1);
   dy = dy(:,1:s2-1);
end;

%======================================================================================
% SAVE DISPLACEMENTS

xv  	= 	mat2vertcol(x);
yv  	= 	mat2vertcol(y);
uxv 	= 	mat2vertcol(dx);
uyv 	= 	mat2vertcol(dy);

folder_for_disp = strcat('/Users/Rose/Dropbox/Postdoc/Algorithms/Singlecell/displacement/',cellfolder,'/'); 

      
 f3 = fopen([folder_for_disp name2 '.dat'],'w');

for i = 1:size(uxv,1),
   fprintf(f3,'%15.5e%15.5e%15.5e%15.5e\n',xv(i), yv(i), uxv(i), uyv(i));   

%     m(i,:)           = [xv(i),yv(i),uxv(i),uyv(i)];
  
end;

%* if floor(str2num(strsec)/10000) < 1
%*             timestr = strcat('0',num2str(strsec));
%*         else
%*             timestr = strcat(num2str(strsec));
%*         end;
        
%  header = '';
%  colnames = '';
%  xlsfilename = strcat(name3,'.xls');
% %  cd(strcat(folder_for_disp,'\Displacement\',cellfolder,'\'))
%       xlswrite(m,header,colnames,xlsfilename);     
% %  cd(strcat(folder_for_disp,'\croppeddata\',cellfolder,'\'))
% 
 fclose(f3);

%======================================================================================

%disp(' '); disp(['  ' num2str(etime(clock,time1),'%10.1f') ' sec, or '...
%      num2str(etime(clock,time1)/60,'%10.1f') ' min.']); disp(' ');

%======================================================================================

