%%
% Settings
TimeLapse = true;

%%
Folder = '/Users/Rose/Dropbox/shared with pooyam/Science/TractionForceMicroscopy';
Date = '20160314_';
Sample = 'MTAC_5mm';
ResFold = [Folder '/Results/' Date '_' Sample];
gif_name = [Folder '/Results/' Date '_' Sample '_Beads.gif'];
prop.pixelsize = 0.22727; % 0.37879; % size of each pixel in micron
prop.young = 5000;

%% Data information
addr = util.gen_addr([Folder '/rawdata/2016_04_14_MTAC_5mm']);
[read_im_ts, read_im_fin, im_props] = Get_Read_Func(addr);
% read_im_ts(Position, TimeStamp, Channel)
% read_im_fin(Position, Channel)

%% Setup
% if exist(gif_file ,'file'), delete(gif_file); end
%20160301-MCF7_sylgard --> 5, 
%20160307_MCF7_PAA --> 4,5,6,10,11
pos =5; % Position
gif_name = [Date '_' Sample '_Pos' num2str(pos) '.gif'];

%% Utils
% util.multi_im_show(read_im_ts(pos, 0, 1), read_im_fin(pos, 1));

%% Crop
title('Crop')
[~,rect] = imcrop(read_im_ts(pos, 0, 2), []);
close, drawnow

%% Adjust size of crop
% % make the image 32x by 32x
% rect = round(rect);
% rect(3:4) = round(rect(3:4)/32)*32-1;
% rect(3) = rect(4);
% crop = @(x) imcrop(x, rect);
%
% Make the images with 415, 415 width
rect = round(rect);
rect(3) = 511;
rect(4) = 511;
crop = @(x) imcrop(x, rect);
%% Read final image
imFB_fin = crop(read_im_fin(pos, 0));
imBigFB_fin = crop(read_im_fin(pos, 1));
imBF_fin = crop(read_im_fin(pos, 2));

%% Run through each time frame
for ts = 0:im_props.max_ts
    % Find and compensate drift
    imFB_org = read_im_ts(pos, ts, 0);
    imFB = crop(imFB_org);
    
    [drift_x, drift_y] = disp_on_blocks(imFB_fin, imFB, round(size(imFB_fin), 1), 0); % TODO size
    imFB = imcrop(imFB_org, rect+[drift_x drift_y 0 0]);
    
    if isempty(imFB) || ~all(size(imFB) == size(imFB_fin))
        error('displacement trim results out of bounds of image');
    end
    
    imBigFB = imcrop(read_im_ts(pos, ts, 1), rect+[drift_y drift_x 0 0]);
    imBF = imcrop(read_im_ts(pos, ts, 2), rect+[drift_y drift_x 0 0]);
    
    % Correlation
    [x, y, dx, dy] = beads_imcorr_bart(imFB_fin, imFB);
    

    %% Traction Force and plots
    
    [xrub, yrub] = Get_Cell_Boundaries(imBigFB);
    %%
    [rmst_iterative,max_traction, theta0, Trace_moment, Uecm, sumforce, prestress,...
        area_cell] = puppy_code(prop, x, y, dx , dy, xrub, yrub, ResFold);
%%    
    rmst = rmst_iterative * 1e12; %RMST traction (Pa);
    max_rmst = max_traction *1e12;
    orientation = theta0*180/pi; % Orientation of principle tractions
    netmoment = -1 * Trace_moment * 1e12;% Net contractile moment (pNm)
    StrainEnergy = Uecm*1e12; %Total strain energy (pJ)
    MaxCumForce = max(sumforce)*1e-3; % Max cumulative force (nN)
    % Prestress = stnanmean(prestress);% Prestress (Pa)
    Area = area_cell; %Projected area of the cell (\mum^2)

%  posMat = [posMat pos];
%  rmstMat = [rmstMat rmst];
%  StrainEnergyMat = [StrainEnergyMat StrainEnergy];
%  AreaMat = [AreaMat Area];
%     save( 'MDA468Syl.mat','posMat', 'rmstMat', 'StrainEnergyMat', 'AreaMat');
    %% GIF animation
% % %     imagesc(sqrt(dx.^2+dy.^2)); axis equal
figure
clear imGIF;
imGIF(:,:,1) = imFB(2:end-2,2:end-2);
imGIF(:,:,2) = imFB_fin(2:end-2,2:end-2);
for i = 1:2
    imshow(imGIF(:,:,i),[])
    if ~isempty(gif_name)
        util.write_gif(gif_name, .5);
    end
end
end

%% Playback GIF animation
% if ispc && ~isempty(gif_file) && exist(gif_file, 'file')
%     winopen(gif_file)
% end
