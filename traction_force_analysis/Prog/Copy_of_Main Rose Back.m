clear all
%% Add path to the program

% read lif images;
[filename, pathname]=uigetfile({'*.lif','Leica Image Format (*.lif)'});
filename=[pathname filename];

for pos = 1:1
    
    im_ts=ci_loadLif(filename,pos);
    imFB_org_ts = im_ts.Image{2};
    imBigFB_ts = im_ts.Image{1};
    imBF_ts = im_ts.Image{3};
    
    im_fin=ci_loadLif(filename,40+pos);
    imFB_fin_ts = im_fin.Image{2};
    imBigFB_fin_ts = im_fin.Image{1};
    imBF_fin_ts = im_fin.Image{3};
    
    %%
    % Settings
    Date = '20160314';
    Sample = 'MTAC_Shear_TF';
    ResFolder = '/Users/Rose/Dropbox/shared with pooyam/Science/TractionForceMicroscopy/Results';
    
    if ~isequal(exist([ResFolder '/Results/' Date], 'dir'),7)
    mkdir(ResFolder,Date)
    end
    prop.pixelsize = 0.22727; % 0.37879; % size of each pixel in micron
    prop.young = 5000;
    %% Data information
    %addr = util.gen_addr([ResFolder '/rawdata/2016_04_14_MTAC_5mm']);
    %[read_im_ts, read_im_fin, im_props] = Get_Read_Func(addr);
    % Example of Get_Read_Func outputs which are also functions:read_im_ts and
    % read_im_fin
    % read_im_ts(Position, TimeStamp, Channel)
    % read_im_fin(Position, Channel)
    %% Utils
    % util.multi_im_show(read_im_ts(pos, 0, 1), read_im_fin(pos, 1));
    %% Crop
    title('Crop')
    [~,rect] = imcrop(imBigFB_ts(:,:,:,1), []);
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
    rect(3) = 607;
    rect(4) = 607;
    crop = @(x) imcrop(x, rect);
    %% Read final image
    %     imFB_fin = crop(read_im_fin(pos, 0));
    %     imBigFB_fin = crop(read_im_fin(pos, 1));
    %     imBF_fin = crop(read_im_fin(pos, 2));
    %% Crop final images
    imFB_fin = crop(imFB_fin_ts);
    imBigFB_fin = crop(imBigFB_fin_ts);
    imBF_fin = crop(imBF_fin_ts);
    
    % get cell boundaries for one t series images
    imBigFB = imcrop(imBigFB_ts(:,:,:,1), rect);
    [xrub, yrub] = Get_Cell_Boundaries(imBigFB);
    %% Run through each time frame
    for ts = 1:size(imFB_org_ts,4)
        % Find and compensate drift
        %imFB_org = read_im_ts(pos, ts, 0);
        
        imFB = crop(imFB_org_ts(:,:,:,ts));
        
        [drift_x, drift_y] = disp_on_blocks(imFB_fin, imFB, round(size(imFB_fin), 1), 0); % TODO size
        imFB = imcrop(imFB_org_ts(:,:,:,ts), rect+[drift_x drift_y 0 0]);
        
        if isempty(imFB) || ~all(size(imFB) == size(imFB_fin))
            error('displacement trim results out of bounds of image');
        end
        
        imBigFB = imcrop(imBigFB_ts(:,:,:,ts), rect+[drift_y drift_x 0 0]);
        imBF = imcrop(imBF_ts(:,:,:,ts), rect+[drift_y drift_x 0 0]);
        
        %% Correlation
        [x, y, dx, dy] = beads_imcorr_bart(imFB_fin, imFB);
        
        %% get cell boundaries at each step  
        %[xrub, yrub] = Get_Cell_Boundaries(imBigFB);
       %% Traction Force and plots
        [rmst_iterative,max_traction, theta0, Trace_moment, Uecm, sumforce, prestress,...
            area_cell] = puppy_code(prop, x, y, dx , dy, xrub, yrub, [ResFolder '/' Date '/' Sample '_position' num2str(pos) '_t' num2str(ts)]);
        %%
        rmst(ts+1) = rmst_iterative * 1e12; %RMST traction (Pa);
        max_rmst(ts+1) = max_traction *1e12;
        orientation(ts+1) = theta0*180/pi; % Orientation of principle tractions
        netmoment(ts+1) = -1 * Trace_moment * 1e12;% Net contractile moment (pNm)
        StrainEnergy(ts+1) = Uecm*1e12; %Total strain energy (pJ)
        MaxCumForce(ts+1) = max(sumforce)*1e-3; % Max cumulative force (nN)
        % Prestress = stnanmean(prestress);% Prestress (Pa)
        Area(ts+1) = area_cell; %Projected area of the cell (\mum^2)
        
        %% Beads GIF animation
%         fig1 = figure;
%         imGIF = zeros(size(imFB(2:end-2,2:end-2),1), size(imFB(2:end-2,2:end-2),2),2);
%         imGIF(:,:,1) = imFB(2:end-2,2:end-2);
%         imGIF(:,:,2) = imFB_fin(2:end-2,2:end-2);
%         gif_name = [ResFolder '/' Date '/' Sample '_position' num2str(pos) '_t' num2str(ts) '_beads.gif'];
%         for i = 1:2
%             imshow(imGIF(:,:,i),[])
%             if ~isempty(gif_name)
%                 util.write_gif(gif_name, .5);
%             end
%         end
%         close(fig1)
        %% GIF animation

        imbead (:,:,:,2*ts-1) = imFB(:,:,1);
        imbead (:,:,:,2*ts) = imFB_fin;
        imFlo (:,:,:,ts) = imBigFB;
        imBright(:,:,:,ts) = imBF;
        %         imutraction (:,:,:,ts) =
        %         imContrac (:,:,:,ts) =
        %         imDisp (:,:,:,ts) =
    end
    
end
%%
for i = 1:size(imbead,4)
    imshow(imbead(:,:,:,i),[]);
    util.write_gif([ResFolder '/' Date '/' Sample '_position' num2str(pos) '_BeadSeries.gif'], 0.05);
end


for i = 1:size(imFlo,4)
    imshow(imFlo(:,:,:,i),[]);
    util.write_gif([ResFolder '/' Date '/' Sample '_position' num2str(pos) '_FluoSeries.gif'], 0.1);
end

for i = 1:size(imBright,4)
    imshow(imBright(:,:,:,i),[]);
    util.write_gif([ResFolder '/' Date '/' Sample '_position' num2str(pos)  '_BrightfieldSeries.gif'], 0.1);
end

%% Playback GIF animation
% if ispc && ~isempty(gif_file) && exist(gif_file, 'file')
%     winopen(gif_file)
% end
