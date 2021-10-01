clc
clear
%% Settings
boundary_drawing_step_iteration = 2;
crop_size = 29 * 32 - 1;

prop.pixelsize = 0.18; %0.22727; % 0.37879; %0.455 %0.18 for 63 and zoom1 size of each pixel in micron (change if you use different objective)

prop.young = 2000;

% Positions which need modifying rect
modify_roi_pos = [];

%% Basic Functions Setup
lif_load = @ci_loadLif.ci_loadLif;

%export_fig = @export_fig.export_fig;

% read lif images;
load('path.mat', 'filename')
[f, pathname] = uigetfile({'*.lif','Leica Image Format (*.lif)'}, 'Select LIF file', filename);
if ~f, warning('No file selected.'); return; end
filename = [pathname f];
save('path.mat', 'filename');

warning('off', 'Images:initSize:adjustingMag');

%% Inits

[pathstr,name,ext] = fileparts(filename);
ResFolder = util.gen_addr(['../Results/' name]);

if ~exist(util.gen_addr(ResFolder), 'dir')
    mkdir(ResFolder)
end

if exist([ResFolder '/ROI.mat'])
load([ResFolder '/ROI.mat'], 'rect', 'xrub', 'yrub');
else
    rect = [];
    xrub = {};
    yrub = {};
end
%% Main loop

for task = 1:2 % task = 1 means user input selections, task = 2 means processing
    
    if task == 1
        fprintf('Getting user crop:\n');
    else
        fprintf('Processing:\n');
    end
        for pos = 48

        if pos>=1 && pos<=70
        prop.young = 2000;
        elseif pos>=73 && pos<=80
        prop.young = 21000;
        else
            disp('error');
        end
        
        get_roi_this_position = (size(rect, 1) < pos || any(modify_roi_pos == pos));
        if task == 1 && ~get_roi_this_position
            continue;
        end
        
        fprintf('\tReading position %d ... ', pos);
        
        [im_ts, npos] = lif_load(filename,66);   %bead plane
        imFB_org_ts = im_ts.Image{1};
        [im_ts, npos] = lif_load(filename,64);
        imBigFB_ts = im_ts.Image{1};
        imBF_ts = im_ts.Image{1};
        
%         im_fin = lif_load(filename, npos/2 + pos);
im_fin = lif_load(filename, 212);
        imFB_fin_org = im_fin.Image{1}(:,:,1);
        %imBigFB_fin_org =im_fin.Image{3}(:,:,1); %im_fin.Image{1};
        %imBF_fin_org = im_fin.Image{3}(:,:,1);
        
        fprintf('Done!\n\t\t');
        
        %% Create a new folder for position
        Pos_folder = util.gen_addr([ResFolder '/Position' num2str(pos) '/']);
        base_gif_name = @(name) [Pos_folder name '.gif'];
        if task == 1
            if ~exist(util.gen_addr(Pos_folder), 'dir')
                mkdir(Pos_folder);
            end
            delete_gif_names = {'BeadSeries', 'FluoSeries', 'BrightfieldSeries'};
            for k = 1:length(delete_gif_names)
                if exist(base_gif_name(delete_gif_names{k}), 'file')
                    delete(base_gif_name(delete_gif_names{k}));
                end
            end
        end
        
        %% Utils
        % util.multi_im_show(read_im_ts(pos, 0, 1), read_im_fin(pos, 1));
        
        if task == 1 && get_roi_this_position
            contin = true;
            while (contin)
                %% Get Crop Coordinates
                title('Crop');
                
                imshow(imBigFB_ts(:,:,:,1));
                h = imrect(gca, [0 0 crop_size crop_size]);
                setResizable(h, false);
                rct = wait(h);
                if isempty(rct), warning('No ROI selected.'); return; end
                close, drawnow;
                rect(pos,:) = round(rct);
                
                [drift_x, drift_y] = disp_on_blocks(imFB_fin_org, imFB_org_ts(:,:,:,1), round(size(imFB_fin_org), 1), 0);
                crop_and_de_drift = @(im) imcrop(im, rect(pos,:) + [drift_x drift_y 0 0]);
                imFB = crop_and_de_drift(imFB_org_ts(:,:,:,1));
                
%                 imFB = imcrop(imFB_org_ts(:,:,:,1), rect(pos,:));
                imFB_fin = imcrop(imFB_fin_org, rect(pos,:));
                % Check de-drift validity
                if isempty(imFB) || ~all(size(imFB) == size(imFB_fin))
                    warning('displacement trim results out of bounds of image');
                else
                    contin = false;
            end
            save([ResFolder '/ROI.mat'], 'rect', 'xrub', 'yrub');
                end
        end
        
        %% Crop final images
        imFB_fin = imcrop(imFB_fin_org, rect(pos,:));
        %imBigFB_fin = imcrop(imBigFB_fin_org, rect(pos,:));
        %imBF_fin = imcrop(imBF_fin_org, rect(pos,:));
        
        %% Run through each time frame
        
        num_ts = size(imFB_org_ts,4);
        
        for ts = [1:num_ts]
            if task == 1 && mod(ts, boundary_drawing_step_iteration) ~= 1
                continue;
            end
            
            if task == 2
                fprintf('%d, ', ts);
            end
            
            % Find and compensate drift
            [drift_x, drift_y] = disp_on_blocks(imFB_fin_org, imFB_org_ts(:,:,:,ts), round(size(imFB_fin_org), 1), 0); % TODO size
            crop_and_de_drift = @(im) imcrop(im, rect(pos,:) + [drift_x drift_y 0 0]);
            
            % get cell boundaries for one t series images
            imBigFB = crop_and_de_drift(imBigFB_ts(:,:,:,ts));
            imBF = crop_and_de_drift(imBF_ts(:,:,:,ts));
            imFB = crop_and_de_drift(imFB_org_ts(:,:,:,ts));
            
            % Check de-drift validity
            if isempty(imFB) || ~all(size(imFB) == size(imFB_fin))
                warning('drift compensation out of bounds for ts = %d', ts);
                continue;
            end
            
            if task == 1 % implicit: mod(ts, boundary_drawing_step_iteration) == 1
                % Get cell boundaries at every boundary_drawing_step_iteration
                [xtemp, ytemp] = Get_Cell_Boundaries(imBigFB, imBF);
                for k = 0:boundary_drawing_step_iteration-1
                    xrub{pos, ts+k} = xtemp;
                    yrub{pos, ts+k} = ytemp;
                end
                save([ResFolder '/ROI.mat'], 'rect', 'xrub', 'yrub');
                continue;
            end
            
            %% Correlation
            [x, y, dx, dy] = beads_imcorr_bart(imFB_fin, imFB);
            
            %% Traction Force and plots
            [rmst_iterative,max_traction, theta0, Trace_moment, Uecm, sumforce, prestress,...
                area_cell, mean_displacement] = puppy_code(prop, x, y, dx , dy, xrub{pos,ts}, yrub{pos,ts}, Pos_folder, ts);
            
            %% Store Results
            rmst(pos,ts) = rmst_iterative * 1e12; %RMST traction (Pa);
            max_rmst(pos,ts) = max_traction *1e12;
            orientation(pos,ts) = theta0*180/pi; % Orientation of principle tractions
            netmoment(pos,ts) = -1 * Trace_moment * 1e12;% Net contractile moment (pNm)
            StrainEnergy(pos,ts) = Uecm*1e12; %Total strain energy (pJ)
            MaxCumForce(pos,ts) = max(sumforce)*1e-3; % Max cumulative force (nN)
            % Prestress = stnanmean(prestress);% Prestress (Pa)
            Area(pos,ts) = area_cell; %Projected area of the cell (\mum^2)
            Displacement (pos,ts)= mean_displacement;
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
            
            util.write_gif(base_gif_name('BeadSeries'), 0.5 + (ts == 1)*1.5, imFB);
            util.write_gif(base_gif_name('BeadSeries'), 0.5, imFB_fin);
            
            util.write_gif(base_gif_name('FluoSeries'), 0.5 + (ts == 1)*1.5, imBigFB);
            util.write_gif(base_gif_name('BrightfieldSeries'), 0.5 + (ts == 1)*1.5, imBF);
            
        end
        fprintf('\n');
        end
    %save([ResFolder '.mat'], 'rmst', 'max_rmst', 'orientation', 'netmoment', 'StrainEnergy', 'MaxCumForce', 'Area','Displacement');

end

save([ResFolder '.mat'], 'rmst', 'max_rmst', 'orientation', 'netmoment', 'StrainEnergy', 'MaxCumForce', 'Area','Displacement');

%% Playback GIF animation
% if ispc && ~isempty(gif_file) && exist(gif_file, 'file')
%     winopen(gif_file)
% end
