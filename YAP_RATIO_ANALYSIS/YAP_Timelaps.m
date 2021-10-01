%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name of file: YAP_Newsha_Timelaps
% Version 2.0

% Purpose: to calculate the ratio of intensity of YAP in nucleus to
% cytoplasm.

% Input variables: 
%   Pos: the position(s) in the lif file to analayze
%   time_interval: time interval between slides in mins
%   crop_size: size of the square ROI to analyze

% created: Nov 29, 2016 13:40 PM
% Author: Rosa Kaviani
% Revisions: In this version, Nucleus segmentation will be performed by
% simple thresholding rather than K-mean clustering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
%% Settings
time_interval = 19; %time interval between slides in mins
crop_size = 800; % size of square ROI to analyze

%% Initialization
lif_load = @ci_loadLif.ci_loadLif;
[name, pathname] = uigetfile({'*.lif','Leica Image Format (*.lif)'}, 'Select LIF file');
if ~name, warning('No file selected.'); return; end
filename = [pathname name];

ResFolder = util.gen_addr([ pwd '/Results/' name(1:end-4)]);

if ~exist(util.gen_addr(ResFolder), 'dir')
    mkdir(ResFolder)
end

for npos=6  % the position(s) to analayze

[im_ts, posi] = lif_load(filename,npos);
imcyto = im_ts.Image{2};
imnuc = im_ts.Image{4};

imcell(:,:,1)=imnuc(:,:,:,1);
imcell(:,:,2)=imcyto(:,:,:,1);
imcell(:,:,3)=0;

imshow(imcell(:,:,:,1),[]);
h = imrect(gca, [0 0 crop_size crop_size]);
setResizable(h, false);
rct = wait(h);
if isempty(rct), warning('No ROI selected.'); return; end
close, drawnow;
rect(npos,:) = round(rct);

%%
fig1 = figure;
hold on
%t=1:size(imcyto,4)
for t=1:size(imcyto,4)
    
    imnuc_crop(:,:,t) = imcrop(imnuc (:,:,1,t), rect(npos,:));
    
    imcyto_crop(:,:,t)= imcrop(imcyto (:,:,1,t), rect(npos,:));
    
    NucMask(:,:,t) = util.MaskNucl(imnuc_crop (:,:,t));
    CytoMask(:,:,t) = util.MaskSegment (imcyto_crop(:,:,t));
    temp = CytoMask(:,:,t);
    temp (NucMask(:,:,t))=0;
    CytoMask(:,:,t) = temp;
    
    
    tempnuc = imcyto_crop(:,:,t);
    tempnuc(~NucMask(:,:,t))=0;
    im_nuc_YAP = tempnuc;
    tempcyto = imcyto_crop(:,:,t);
    tempcyto(~CytoMask(:,:,t))=0;
    im_cyto_YAP = tempcyto;
    
    YAP_intensity_nuc(:,:,t) = sum(sum(nonzeros(im_nuc_YAP)));
    YAP_size_pixels(:,:,t) = length(nonzeros(im_nuc_YAP));
    YAP_intensity_cyto(:,:,t) = sum(sum(nonzeros(im_cyto_YAP)));

    Nuc2Cyto (npos,t) = (YAP_intensity_nuc(:,:,t)./size(nonzeros(im_nuc_YAP),1))/(YAP_intensity_cyto(:,:,t)...
        /size(nonzeros(im_cyto_YAP),1));
    NucYAP(npos, t) = YAP_intensity_nuc(:,:,t);
    NucYAP_size(npos,t) = YAP_size_pixels(:,:,t);
    CytoYAP(npos, t) = YAP_intensity_cyto(:,:,t);
    
    %Area plots 
    
    % plot
    time = (0:size(Nuc2Cyto(npos,:),2)-1).*time_interval;
    s1 = subplot(2,2,[1; 3]);
    %set(s1, 'pos', p1);
    h1 = plot(time,Nuc2Cyto(npos,:),'bo','MarkerFaceColor','b','markers',3);
    xlabel('time (min)');
    ylabel ('YAP intensity nucleus/cytoplasm')
    %h1 = plot(1:size(NucYAP,2),NucYAP,'ro','MarkerFaceColor','r');
   % xlim([0 (size(imcyto,4)-1)*time_interval]);
   xlim([0 t*time_interval]);
    %ylim([0 0.5*10e4]);
    ylim([0 3]);
    s2 = subplot(2,2,2);
     p2 = get(s2, 'pos'); %[left, bottom, width, height] 
    p2(3) = p2(3)*1.10;
    p2(4) = p2(4)*1.10;
    set(s2, 'pos', p2);
    imshow(im_cyto_YAP,[]);
    s3 = subplot(2,2,4);
      p3 = get(s3, 'pos'); %[left, bottom, width, height] 
    p3(3) = p3(3)*1.10;
    p3(4) = p3(4)*1.10;
    set(s3, 'pos', p3);
    imshow(im_nuc_YAP,[]);
    pause(0.5);
    frame = getframe(1);
    clear tempnuc
    clear tempcyto
        
    
   util.write_gif([ResFolder '\Nuc2CytoRealtime.gif'], 0.5 + (t == 1)*1.5, frame.cdata);

end

 time = (0:size(NucYAP_size(npos,:),2)-1).*time_interval;
    NucYAP_size(npos,:) = 0.227^2*NucYAP_size(npos,:);
    f5 = figure; 
    figure(f5)
    hold on 
    plot(time, NucYAP_size(npos,:),'bo','MarkerFaceColor','b','markers',3);
    ylim([0 500]);
   xlabel('time (min)');
    ylabel ('Nucleussize');
    hold off 
   
  
   [R P]= corrcoef(Nuc2Cyto(npos,:),NucYAP_size(npos,:)) 
  
    f6 = figure; 
    figure(f6)
    hold on
 % plot(time,R);

    

end


%save([ResFolder 'YAP2CYTO.mat'], 'Nuc2Cyto','NucYAP','CytoYAP','R','P');

save([ResFolder 'YAP2CYTO.mat'], 'Nuc2Cyto','NucYAP','CytoYAP');


