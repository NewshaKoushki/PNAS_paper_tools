%******************************************************************** 
% Functions:
%    Detect mitochondria movement without twisting field;
%    Used by Ramaswamy Krishnan to separate disp synchronous with the applied
%    magetic field
% 
% Inputs:
%    Time period for Mitochondria data selected; Gauss: 0 or not
%   
% Outputs:
%    1. Displacement figures and movies;
%    2. Consective displacements after excluding drifting and synchronous movement;
%    3. Mean Square Displacement (MSD) with time course;
%
% Displacement detection algorithm modified by Stanley Hu on 07/21/2005
% Adapted by Ramaswamy Krishnan 09/21/05
%********************************************************************

clear all; close all;
no_of_cycles     = 1;
images_per_cycle = 10;
base=[];
%std_cutoff      = 0.002;
std_cutoff_rate  = 0.002;
snr_cutoff       = 1;
gray_threshold   = 0.35;%0.45; %threshold for identifying image including only beads;
weight_threshold = 3/10; %threshold for identifying a true bead;
bead_width       = 70;%24;%40;%36; %bead diameter in pixel,even number;
bead_weight      = 520;%130;%520; %30x30=900 pixels;
amp_threshold    = 0; %0.04; %pixel;

time_digit       = 5;
current_digit    = 4;

new_machine      = 1;
frequency_hz     = 0.75; %frequency of sine waveform, Hz;
um_per_pixel     = 0.3263; %for 20x;
exposure_time    = 0.1;   %sec.

 folds = {'Cell29'};
% folds = {'Cell17'};

for fi1 = 1 : length(folds),
    cellfolder	= folds{fi1}
    eval(['Syn_mvmt']); 

%Selected window size, evaluated from Syn_mvmt;
col_start = 1;
col_end   = j2-j1+1;
row_start = 1;%20;11
row_end   = i2-i1+1;
no_of_points = 1; % traction is evaluated every 'no_of_points'

%----------------------------------------------------------------------------------------------------------
image_directory = strcat('C:\C2C12\croppeddata\Syncronized\',cellfolder);
% prog_directory= strcat('C:\C2C12\croppeddata\',cellfolder);
save_directory = strcat('C:\C2C12\Displacement\Syncronized\',cellfolder);
hand_select = 1; %flag for hand-selecting an interested area;
first_cell = 1; 
last_cell = length(min_time); 

% 1______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1_______________________________1___________________________1____________________________1____________________________1____________________________1____________________________1____________________________1____________________________1____________________________1____________________________1____________________________1____________________________1____________________________1____________________________1____________________________1____________________________1____________________________1____________________________1____________________________1____________________________1____________________________1____________________________1____________________________1____________________________1____________________________1____________________________1____________________________1____________________________1____________________________1____________________________1____________________________1____________________________1____________________________1____________________________1____________________________1____________________________1____________________________1____________________________1____________________________1____________________________1____________________________1____________________________1____________________________1____________________________1____________________________1____________________________1____________________________1____________________________1____________________________1____________________________1____________________________1____________________________1____________________________1____________________________1____________________________1____________________________1____________________________1____________________________1____________________________1____________________________1____________________________1____________________________1_______________________________1____________________________ 1____________________________ 1____________________________ 1____________________________ 1____________________________ 1____________________________ 1____________________________ 1____________________________ 1____________________________1____________________________ 1____________________________ 1____________________________ 1____________________________ 1____________________________ 1____________________________ 1____________________________ 1____________________________ 1____________________________ 1________________________11__1________________________11__1________________________11__1________________________11__1________________________11__1________________________11__1________________________11__1________________________11__1________________________11__1________________________11__1________________________11__1________________________11__1________________________11__1________________________11__1________________________11__1________________________11__1________________________11__1________________________11__1________________________11__1________________________11__1________________________11__1________________________11__1________________________11__1________________________11__1________________________11__1________________________11__1________________________11__1________________________11__1________________________11__1________________________11__1________________________11__1________________________11__1________________________11__1________________________11__1________________________11__1________________________11__1________________________11__1________________________11__1________________________11__1________________________11__1________________________11__1________________________11__1________________________11__1________________________11__1________________________11__1________________________11__1________________________11__1________________________11__1________________________11__1________________________11__1________________________11__1________________________11__1________________________11__1________________________11__1________________________11__1________________________11__1________________________11__1________________________11__1________________________11__        
 base =     ['cropCell1'; 'cropCell1'; 'cropCell1'; 'cropCell1'; 'cropCell1'];
%Cell 2
% min_time = [ 07985; 08738; 12115; 17633; 20907];
% max_time = [ 07999; 11910; 15057; 20573; 24617];       
% Cell 4
% min_time = [ 07911; 08682; 11812; 15207; 18408];
% max_time = [ 07972; 11792; 14752; 18205; 21886];                          
%Cell 10
% min_time = [ 02695; 03359; 06407; 09569; 12648];
% max_time = [ 02717; 06299; 09397; 12565; 15654];       
%Cell 15
% min_time = [ 03773; 04163; 07378; 10567; 13703];
% max_time = [ 03780; 07165; 10374; 13673; 16603];       



%----------------------------------------------------------------------------------------------------------



%check time periods
for j = 1:length(min_time),
    if (min_time(j)>=max_time(j)) ,
        disp('Time period error! j=');
        j
        return;
    end;
end;
    

base = cellstr(base);

shift_step = 0.05 ; %pixel
shift_max  = 4/shift_step +1;
shift_center = 4/shift_step/2+1;

% compute matrix "shift" for shift in Fourier space
% template  8 for a resolution of 8 pixels. template = 5 for a 3 pixel resolution
template = 5; %size of template in pixels = (2*template)+1
% margin = 5 for a resolution of 8 pixels. Margin = 3 for a 3 pixel resolution
margin   = 3; %paramater for transformation from image to analysis coordinate system
T = 2*template+5;
Nyquist = (T/2+1);    
disp('if shift variable does not exist, save shift variable into a Matlab file');
save_file = strcat(save_directory,'\','shift_finer_',num2str(template,'%03d'),'_Step',num2str(shift_step*1000,'%03d'),'.mat');
%save_file = strcat(prog_directory,'shift_finer_',num2str(template,'%03d'),'.mat');
if exist(save_file),
    load (save_file);
else
	for axi = 1:shift_max,
        axi/shift_max*100
        for ayi = 1:shift_max,
            ax = (axi-shift_center)*shift_step; %shift in x and y direction is maximally 2 pixels in steps of 1/10th of a pixel
            ay = (ayi-shift_center)*shift_step;
			for x = 1:T,
				for y = 1:T,                
                    if (x<=Nyquist) & (y<=Nyquist),
                        shift(ayi,axi,y,x) = exp(-i*2*pi/T*(ax*(x-1)+ay*(y-1)));
                    elseif (x>Nyquist) & (y<=Nyquist),
                        shift(ayi,axi,y,x) = exp(-i*2*pi/T*(-ax*(T-(x-1))+ay*(y-1)));
                    elseif (x<=Nyquist) & (y>Nyquist),
                        shift(ayi,axi,y,x) = exp(-i*2*pi/T*(-ay*(T-(y-1))+ax*(x-1))); 
                    else
                        shift(ayi,axi,y,x) = exp(-i*2*pi/T*(-ax*(T-(x-1))-ay*(T-(y-1))));     
                    end;
                end;            
            end;
        end;
	end;
    save (save_file, 'shift');
end;    

clear bw2;

for cell_no = first_cell:last_cell,
      
    close all;
    clear b b1 b2 bw1 bw3 b_bead b_movie covar snr current currentX currentY currentZ bw coord time sinus cosinus final_x final_y model_x model_y 
    clear F shift_circle shift_bound template_circle template_bound
	%find all bitmap files in current directory
    directory = char(base(cell_no));   
    k = findstr(directory,'_p');
    if ~isempty(k),
        directory = directory(1:k-1);	    
    end;
%    directory = strcat(image_directory,directory,'\');
    directory = strcat(image_directory,'\');
    cdirectory = strcat('cd (',char(39), directory, char(39),')');     
	eval(cdirectory); 
	file_list=dir;
	file_nr = [];
	for m = 1:length(file_list)
        k = findstr(file_list(m).name,'.tif');
        if (k>10),
            file_nr = [file_nr m];
        end;
	end;
	file_array = char(file_list(file_nr).name);
	clear file_list file_nr;
	
	%sort file_array
	time=[];
	for m=1:length(file_array(:,1)),
        k = findstr(file_array(m,:),'_');
        time_string = file_array(m,k+1:k+time_digit);
        time = [time str2num(time_string)];  
	end;
	[y,k]=sort(time);
	file_array = file_array(k,:);

    % find first  bright image;
	time     = [];
	current  = [];
%     currentX = [];
%     currentY = [];
%     currentZ = [];
%     
	for m=1:length(file_array(:,1)),
        k = findstr(file_array(m,:),'_');
        time_string = file_array(m,k+1:k+time_digit);
        time = [time str2num(time_string)];
        current_string = file_array(m,k+time_digit+1:k+time_digit+1+current_digit);
%         currentX       = [currentX str2num(current_string)];    
%         current_string = file_array(m,k+time_digit+current_digit+2:k+time_digit+2*(current_digit+1));
%         currentY       = [currentY str2num(current_string)];
%         current_string = file_array(m,k+time_digit+2*(current_digit+1)+1:k+time_digit+3*(current_digit+1));
%         currentZ       = [currentZ str2num(current_string)];      
        current        = [current str2num(current_string)];    
	end;
    
    frequency_hz = freq_s(cell_no);
   
	
    %find the first Mitochondria image;
    if (gauss_s(cell_no) <= 0),
        m_min = min(find(time>=min_time(cell_no)));  
		m_max = max(find(time<=max_time(cell_no)));
    else
        m_min = min(find(time>=min_time(cell_no)))-1;
%        m_min = min(find(time>=min_time(cell_no)));
        m_max = max(find(time<=max_time(cell_no)));
        time_limit = time(m_min:m_max);
        current_limit = current(m_min:m_max);
        cross_points_index = [];
        for index = 2:length(time_limit)-1,
            if ((current_limit(index) > 0) & (current_limit(index-1) < 0) & (current_limit(index) < current_limit(index+1))),
                cross_points_index = [cross_points_index index];
            end;           
        end; %for
              
        if (length(cross_points_index) >= 2),   
            b_begin_index = cross_points_index(1);        
            for index = 2:length(cross_points_index),
                if (mod((cross_points_index(index)-b_begin_index), images_per_cycle)==0),
                    m_min_b = b_begin_index;
                    m_max_b = cross_points_index(index);        
                else
                    b_begin_index = cross_points_index(index);
                end;
             end; %for
        end;
        
        min_time_b = time_limit(m_min_b);
        max_time_b = time_limit(m_max_b);      
        m_min = min(find(time>=min_time_b));
        m_max = max(find(time<=max_time_b));
        m_max = m_min + 10*floor((m_max-m_min)/10)-1    
        no_of_cycles = (m_max+1-m_min)/10;
      
		if (m_max+1-m_min)<images_per_cycle,
            disp('Error: the times for the beginning and end that you have specified');
            disp('are less than the number of images per twisting cycle!');
            return;
		end;
            
        
        
    end;
          
    %make a movie for images;
    figure(2); 
    rect = [100, 50, 300, 300];
    fig=figure(2);
    set(fig,'NextPlot','replace','Position',rect,'Visible','off');
%     movie_file = strcat(save_directory,'movie_',char(base(cell_no)),'_',num2str(time(m_min)),'_',num2str(time(m_max)),'.avi');
%     mov = avifile(movie_file,'Compression','none','FPS',10);
    
    
    %make a movie for displacements;
    figure(4); 
    rect = [600, 50, 300, 300];
    fig=figure(4);
    set(fig,'NextPlot','replace','Position',rect,'Visible','off');
%     movie_file = strcat(save_directory,'movie_Xdisp_',char(base(cell_no)),'_',num2str(time(m_min)),'_',num2str(time(m_max)),'.avi');
%     movX = avifile(movie_file,'Compression','none','FPS',1);

    figure(5); 
    rect = [600, 400, 300, 300];
    fig=figure(5);
    set(fig,'NextPlot','replace','Position',rect,'Visible','off');
%     movie_file = strcat(save_directory,'movie_Ydisp_',char(base(cell_no)),'_',num2str(time(m_min)),'_',num2str(time(m_max)),'.avi');
%     movY = avifile(movie_file,'Compression','none','FPS',1);
    
    figure(6); 
    rect = [100, 400, 300, 300];
    fig=figure(6);
    set(fig,'NextPlot','replace','Position',rect,'Visible','off');
%     movie_file = strcat(save_directory,'movie_TotalDisp_',char(base(cell_no)),'_',num2str(time(m_min)),'_',num2str(time(m_max)),'.avi');
%     movAmp = avifile(movie_file,'Compression','none','FPS',1);
%   
	
   
    
	%********************************************************************************************************************
	m_first = m_min-1;  %-1: To get the displacement data of exact entire cycles (e.g. nx10 points) with non-zero Gauss;
    m_last = m_max;   
	covar = [];
    bw = [];
    b1 = [];
	number = 0  % number of files analyzed
    

    % load and filter first image (b1);
    
    first_file = deblank(file_array(m_first,:))        		
	[b map] = imread(first_file,'tiff');
	b = double(b);%(91:191,41:141));
%     b = b ./ max(max(b));
%     [y,x]=size(b);
%     Nyquistx = (x/2+1);
%     Nyquisty = (y/2+1);
%     bfft2=fft2(b);
%     bfft2(Nyquisty,Nyquistx)=0; %take out highest frequency
%     b=real(ifft2(bfft2));
%     clear bfft2;     
    b1 = b;
    b1 = b1 ./ max(max(b1));
    figure(1); clf; subplot(3,3,1); imagesc(b1); title(first_file); colorbar; drawnow;
    figure(2);imagesc(b1); text(10,10,'0.0','Color','w'); colormap(gray);
%     F=getframe(gca);
%     mov = addframe(mov,F);

    min_b1=min(min(b1));
 
    
    [row,col] = size(b1);
    new_col = floor((col-margin-2-template)/template);
    new_row = floor((row-margin-2-template)/template);
 
    
    covar=zeros(new_row,new_col,2);
    for xx = 1:new_col                        
        for yy = 1:new_row            
            xxx = margin + (xx * template); % center position for x
            yyy = margin + (yy * template); % center position for y 
            ima = b1((yyy-template):(yyy+template),(xxx-template):(xxx+template)); 
            %imb = b_bead((yyy-template):(yyy+template),(xxx-template):(xxx+template));  
            covar(yy,xx,1)=std(std(ima)); 
            covar(yy,xx,2)=mean(mean(ima));
            
        end;
    end; 

        
    if isempty(bw),  
        figure(1); subplot(3,3,3);imagesc(covar(:,:,1),[0 0.02]); title('std'); colorbar; 
        std_cutoff = min(min(covar(:,:,1))) + std_cutoff_rate * (max(max(covar(:,:,1))) - min(min(covar(:,:,1))));  
        %std_cutoff = 0.00212215;
        bw = im2bw(covar(:,:,1)/std_cutoff,1);

        %bw = im2bw(covar(:,:,1),std_cutoff);
        bw = bwmorph(bw,'clean'); 
        bw = bwmorph(bw,'fill');
        bw1 = bwmorph(bw,'bridge');
        %bw1 = bwmorph(bw1,'majority',Inf);
        %bw = bwmorph(bw,'close');
        %bw = bwmorph(bw,'dilate');
        bw = bw | bw1;
        figure(1); subplot(3,3,2);imagesc(bw);title('mask'); colorbar;
                 
        b2=b1;
        d_black=floor(template/2)+1;
        for xx = 1:new_col                        
            for yy = 1:new_row            
                xxx = margin + (xx * template); % center position for x
                yyy = margin + (yy * template); % center position for y                  
                if bw(yy,xx)<1,
                    b2((yyy-d_black):(yyy+d_black),(xxx-d_black):(xxx+d_black))=min_b1;
                end;                    
            end;
        end;            
        figure(1);subplot(3,3,1);imagesc(b2); title(first_file); colorbar; drawnow;
    end;
    
%   ----------------------------------------------------------------------------------------------------------------

%   Specify the x and y coordinates where the displacements will be evaluated
%   The x and y coordinates will be read from the trypsindisp file for a particular cell

%       x_mito_selected = [32 64 96 128 160 192 224];
%       y_mito_selected = [32 64 96 128 160 192 224];
    
     cd(strcat('C:\C2C12\Displacement\',cellfolder)); pwd;
     d2 = dir; si = 1;
     d2b(1)=d2(1);			%d2b is part of d2 and will include [.][..] and original cell files
     d2b(2)=d2(2);
      for d2i = 3 : length(d2),
      if ~isempty(findstr(d2(d2i).name,'trypsindisp')) & ...		% name has to contains 'trypsindisp'
         ~isempty(findstr(d2(d2i).name,'.xls')), 			%name has to contains '.xls'
      		d2b(si+2)	=d2(d2i);
		      name 		= d2(d2i).name
              si = si + 1;
          end;
      end;

      % ~isempty(findstr(d2(d2i).name,num2str(second(cell_no)))) & ...% name identifies first file in freq. loop

      selectpts = xlsread(strcat('C:\C2C12\Displacement\',cellfolder,'\',name)); 
     
      si =0;
      for row = 1: size(selectpts,1),
          
          % delete the edge points
          if selectpts(row,1) ~= selectpts(1,1) & selectpts(row,1) ~= selectpts(length(selectpts),1) &...
              selectpts(row,2) ~= selectpts(1,1) & selectpts(row,2) ~= selectpts(length(selectpts),1),
              si = si+1;
              x_mito_selected(si) = [selectpts(row,1)];
              y_mito_selected(si) = [selectpts(row,2)];
          end %(if statement)
          end; %(for statement)
%      

%   ----------------------------------------------------------------------------------------------------------------       
       
    if hand_select & (cell_no==first_cell),  %'1' means hand-selecting the interested area;  
        
        figure(3);imshow(b1); set(3,'Name',strcat('      ',char(base(cell_no)), ',Click the interested Mito!' ));
        hold on;
%        [x_mito_selected,y_mito_selected,I_mito] = impixel(b1);
        I_mito = impixel(b1,x_mito_selected,y_mito_selected);
        
        x_mito = (x_mito_selected-margin)/template;
        y_mito = (y_mito_selected-margin)/template;
        bw2    = zeros(new_row,new_col);
        
        for k = 1:length(x_mito),

            bw2(round(y_mito(k)), round(x_mito(k))) = 1;
            plot(x_mito_selected(k),y_mito_selected(k),'r*');
            h = text(x_mito_selected(k),y_mito_selected(k),num2str(k));
            set(h,'Color','g');
        end;
 
        save_file = strcat(save_directory,'\','Selectedpoints_',char(base(cell_no)),'.jpg');
        saveas(3,save_file);
       
        
        figure(31); imshow(bw2);
        
        
%         figure(3);imshow(b1); set(3,'Name',strcat('      ',char(base(cell_no)), ',      Select the interested area!' ));        
%         [bw2_selected, xrub_selected, yrub_selected] = roipoly(b1);   %hand-select the area of interest;
%         xrub = (xrub_selected-3)/template;
%         yrub = (yrub_selected-3)/template;
%         bw2_weight  =zeros(new_row,new_col);
%         for xx = 1:new_col                        
%             for yy = 1:new_row            
%                 xxx = 3 + (xx * template); % center position for x
%                 yyy = 3 + (yy * template); % center position for y                  
%                 ima = bw2_selected((yyy-template):(yyy+template),(xxx-template):(xxx+template)); 
%                 bw2_weight(yy,xx)=sum(sum(ima));               
%             end;
%         end;        
%         bw2 = im2bw(bw2_weight,0+eps);      
%         save_file = strcat(image_directory,'Data\','SelectedArea_bw_',char(base(cell_no)),'.mat');
% 	    save (save_file, 'bw2', 'xrub', 'yrub');        
%         close(3);
        

    end;
   
    
    if sum(sum(bw2)) ~= 0, 
       bw = bw | bw2;
    end;
    
    
    clear time_buffer sum_moveX sum_moveY sum_moveTotal avg_moveX avg_moveY avg_moveTotal
    
    deltaX = [];
    deltaY = [];
    
    sum_moveX     = zeros(new_row,new_col);
    sum_moveY     = zeros(new_row,new_col);
    sum_moveTotal = zeros(new_row,new_col);
    avg_moveX     = [];
    avg_moveY     = [];
    avg_moveTotal = [];
    
    
    % go through the next images
% 	m_next = m_first+1;
 
    time_buffer = [];
    for m_next = (m_first +1):step(cell_no):m_last,    
    
% 	while m_next <= m_last,
        b2 = [];
        number=number+1
        time_buffer = [time_buffer; time(m_next)];
        
        % load and filter next image (b2);
       
            next_file=[];
            next_file = deblank(file_array(m_next,:));
              cd(image_directory);
		    [b_original map] = imread(next_file,'tiff');
            b = b_original(row_start:row_end, col_start:col_end);
		    b = double(b);%(91:191,41:141));
%             b = b ./ max(max(b));
%             [y,x]=size(b);
%             Nyquistx = round(x/2+1);
%             Nyquisty = round(y/2+1);
%             bfft2=fft2(b);
%             bfft2(round(Nyquisty),round(Nyquistx))=0; %take out highest frequency
%             b=real(ifft2(bfft2));
%             clear bfft2;

        b2 = b;
      
        b2 = b2 ./ max(max(b2));
        figure(1);subplot(3,3,2);imagesc(b2); title(next_file);  colorbar; drawnow;   
        figure(2);imagesc(b2); text(10,10,num2str(time(m_next)/10,3),'Color','w'); colormap(gray);
%         F=getframe(gca);
%         mov = addframe(mov,F);

        for xx = 1:new_col                
            set(1,'Name',strcat(num2str(number),'done;   Total:', num2str(m_max-m_min+1)));                     
            for yy = 1:new_row
                %axi=shift_center; ayi=shift_center; %starting position
                if bw(yy,xx,1)>0.9,  
                    xxx = margin + (xx * template); % center position for x
                    yyy = margin + (yy * template); % center position for y 
                    ima = b1(yyy-template-2:yyy+template+2,xxx-template-2:xxx+template+2);
                    imb = b2(yyy-template:yyy+template,xxx-template:xxx+template);      
                    fta = fft2(ima);
                    
                    %Searching Algorithm by Stanley Hu
                    ay_series = []; ax_series = [];
                    imrms(1:shift_max,1:shift_max)=+Inf;
                    axi = shift_center; ayi = shift_center; %starting position
                    while (1),
                         ay_series =[ay_series;ayi];ax_series =[ax_series;axi];

                        
                        if (imrms(ayi,axi) == +Inf),
                            imrms(ayi,axi) = fs2D(T,ayi,axi,shift,fta,ima,imb);
                        end;
                        if (imrms(ayi,axi-1) == +Inf),
                            imrms(ayi,axi-1) = fs2D(T,ayi,axi-1,shift,fta,ima,imb);
                        end;
                        if (imrms(ayi-1,axi) == +Inf),
                            imrms(ayi-1,axi) = fs2D(T,ayi-1,axi,shift,fta,ima,imb);
                        end;
                        if (imrms(ayi,axi+1) == +Inf),
                            imrms(ayi,axi+1) = fs2D(T,ayi,axi+1,shift,fta,ima,imb);
                        end;
                        if (imrms(ayi+1,axi) == +Inf),
                            imrms(ayi+1,axi) = fs2D(T,ayi+1,axi,shift,fta,ima,imb);
                        end;
 
%                         figure(21);subplot(1,2,1);plot(ax_series,ay_series,'*b'); drawnow;
%                         subplot(1,2,2);imagesc(imrms); colorbar;drawnow;
%                         pause;                        
                        
                        local_imrms = [Inf imrms(ayi-1,axi) Inf;
                                       imrms(ayi,axi-1) imrms(ayi,axi) imrms(ayi,axi+1);
                                       Inf imrms(ayi+1,axi) Inf];
                        [local_minY local_minX] = find(local_imrms == min(min(local_imrms)));
                        if ~isempty(find(local_minX.*local_minY==4)), %(2,2)?
                            break; %the match found
                        end;
                        ayi = ayi + local_minY(1) - 2;
                        axi = axi + local_minX(1) - 2;
                        if ((ayi == 1)|(ayi == shift_max)|(axi == 1)|(axi == shift_max)),
                            %axi = shift_center; ayi = shift_center;
                            break;%searching to the edge
                        end;  
                    end; %while
                                           
             
					ay = (ayi-shift_center)*shift_step;
                    ax = (axi-shift_center)*shift_step;
                    coord(number,yy,xx,1) = ax;                            
                    coord(number,yy,xx,2) = ay; 
                                          
                    
                else
                    coord(number,yy,xx,1) = 0; coord(number,yy,xx,2) = 0; 
                    
                end;   

            end;
        end;   
        
        
        clear dd;
        
        %displacement (x) for current image pair
        for yy=1:new_row,
            for xx=1:new_col,
                dd(yy,xx) = coord(number,yy,xx,1);
            end;
        end;
        nonzero = find(dd ~= 0);
        shiftX = mean(mean(dd(nonzero))); 
        sum_moveX(:,:) = sum_moveX(:,:) + dd(:,:);
        avg_moveX = [avg_moveX sum(sum(sum_moveX))/sum(sum(bw))];
        
        figure(1); subplot(3,3,4);  imagesc(dd); title('Step Displacement X'); colorbar; 
        
        figure(4); cla; imagesc(sum_moveX,[-max2(abs(sum_moveX)) max2(abs(sum_moveX))+eps]);  
        title(strcat('Total Displacement X:', num2str(time(m_next)/10,3), 'sec')); colorbar; 
%         F=getframe(gcf);
%         movX = addframe(movX,F);
        
        
        clear dd;
        %displacement (y) for current image pair
        for yy=1:new_row,
            for xx=1:new_col,
                dd(yy,xx) = coord(number,yy,xx,2);
            end;
        end;
        nonzero = find(dd ~= 0);
        shiftY = mean(mean(dd(nonzero))); 
        sum_moveY(:,:) = sum_moveY(:,:) + dd(:,:);
        avg_moveY      = [avg_moveY sum(sum(sum_moveY))/sum(sum(bw))];
        
        figure(1); subplot(3,3,5);  imagesc(interp2(dd,0,'cubic')); title('Step Displacement Y'); colorbar;
        figure(5); cla; imagesc(sum_moveY,[-max2(abs(sum_moveY)) max2(abs(sum_moveY))]); 
        title(strcat('Total Displacement Y:', num2str(time(m_next)/10,3), 'sec')); colorbar; 
%         F=getframe(gcf);
%         movY = addframe(movY,F); 
                
        
        figure(1); subplot(3,3,7);  imagesc(sum_moveX(:,:)); title('Total Displacement X'); colorbar;
        figure(1); subplot(3,3,8);  imagesc(sum_moveY(:,:)); title('Total Displacement Y'); colorbar;
        
        sum_moveTotal(:,:) = (sum_moveX(:,:).^2.0 + sum_moveY(:,:).^2.0) .^ 0.5;
        avg_moveTotal      = [avg_moveTotal sum(sum(sum_moveTotal))/sum(sum(bw))];
        figure(1); subplot(3,3,9);  imagesc(sum_moveTotal(:,:)); title('Total Displacement-Amp'); colorbar;
         
        
        figure(6); cla; imagesc(sum_moveTotal); 
        title(strcat('Total Displacement-Amp:', num2str(time(m_next)/10,3), 'sec')); colorbar; 
%         F=getframe(gcf);
%         movAmp = addframe(movAmp,F); 
                
        
        if  number >= 1,
           figure(1); subplot(3,3,6); title('Avg Disp (X-red,Y-blue,Total-black)'); hold on;
           plot(time_buffer(:)/10, avg_moveX(:),'ro-');
           plot(time_buffer(:)/10, avg_moveY(:),'bo-'); 
           plot(time_buffer(:)/10, avg_moveTotal(:),'ko-');
           %axis([0.9*time(m_min+1) 1.1*time(m_min+number) -Inf Inf]);
           hold off;
        end;
   
        deltaX_vect = [];
        deltaY_vect = [];
        for k = 1:length(x_mito),
                

            deltaX_vect = [deltaX_vect (coord(number,round(y_mito(k)), round(x_mito(k)),1)-shiftX)*um_per_pixel];
            deltaY_vect = [deltaY_vect (coord(number,round(y_mito(k)), round(x_mito(k)),2)-shiftY)*um_per_pixel];
%            deltaX_vect = [deltaX_vect (coord(number,round(y_mito(k)), round(x_mito(k)),1)-shiftX)*um_per_pixel];
%            deltaY_vect = [deltaY_vect (coord(number,round(y_mito(k)), round(x_mito(k)),2)-shiftY)*um_per_pixel];
        end;
        
        deltaX = [deltaX; deltaX_vect];
        deltaY = [deltaY; deltaY_vect];
       
        
        %m_next = m_next + 1;
        b1=b2; 
        first_file = [];
        first_file = next_file;        
            
	end; %while m_next <= m_last,	
    
%     mov    = close(mov);
%     movX   = close(movX); 
%     movY   = close(movY);
%     movAmp = close(movAmp);

    
    figure(1); drawnow;
    save_file = strcat(save_directory,'\','fig_',char(base(cell_no)),'_',num2str(time(m_min)),'_',num2str(time(m_max)),'.jpg');
    print('-djpeg90', save_file);   
    
    figure(3); subplot(2,2,1);  hold on;
    plot(time_buffer(:)/10, avg_moveX(:),'bo-'); ylabel('Avg Displacement X'); 
    title(strcat('Average movement is:  ',num2str(mean(avg_moveX(:)))));
    hold off;
    figure(3); subplot(2,2,2);  hold on;
    plot(time_buffer(:)/10, avg_moveY(:),'bo-'); ylabel('Avg Displacement Y'); 
    title(strcat('Average movement is:  ',num2str(mean(avg_moveY(:)))));
    hold off;
    figure(3); subplot(2,2,3);  hold on;
    plot(time_buffer(:)/10, avg_moveTotal(:),'bo-'); ylabel('Avg Displacement-Amp'); 
    title(strcat('Average movement is:  ',num2str(mean(avg_moveTotal(:)))));
    hold off;
    
    save_file = strcat(save_directory,'\','fig_',char(base(cell_no)),'_',num2str(time(m_min)),'_',num2str(time(m_max)),'.fig');
    saveas(3,save_file);  
    close(3);
    
    save_file = strcat(save_directory,'\','coord_',char(base(cell_no)),'_',num2str(time(m_min)),'_',num2str(time(m_max)),'.mat');
	save (save_file, 'coord');
 

        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %save the spatial positions and corresponding displacements summed over prescribed 'no_of_points'
    %X Y DELTX DELTY
 
     si =0;
       for row = 1: size(selectpts,1),
            if selectpts(row,1) ~= selectpts(1,1) & selectpts(row,1) ~= selectpts(length(selectpts),1) &...
              selectpts(row,2) ~= selectpts(1,1) & selectpts(row,2) ~= selectpts(length(selectpts),1),
              si = si+1;
               DX_trypsin(si) = [selectpts(row,3)];
               DY_trypsin(si) = [selectpts(row,4)];
           end;   %(if statement)
       end
       
       DX = [];
    DY = [];          
     m = [];
     
	disp('store data as txt-file');
   	data_file = strcat(save_directory,'\DeltaDisplacement_',char(base(cell_no)),'_',num2str(time(m_min)),'_',num2str(time(m_max)),'.csv');
    fid1 = fopen(data_file,'w');
   
 	for i = 1: floor(length(time_buffer)/no_of_points),
                     DX(i,:)            = sum(deltaX(1:(no_of_points*i),:),1);
                     DY(i,:)            = sum(deltaY(1:(no_of_points*i),:),1);
        
        for j = 1:length(x_mito),    
             DispXtxt        = num2str(deltaX(i,j));
             DispYtxt        = num2str(deltaY(i,j));
              m(j,:)         = [x_mito_selected(j),y_mito_selected(j),(deltaX(i,j)),(deltaY(i,j))];
         end;
  
          header = '';
        colnames = '';

%     Add a leading zero to time points with 4 digits (so that windows arranges files correcly by name)
        if floor(time_buffer(i*no_of_points))/10000 < 1
            timestr = strcat('0',num2str(time_buffer(i*no_of_points)));
        else
            timestr = num2str(time_buffer(i*no_of_points));
       end;

      xlsfilename = strcat('trypsindisp_',num2str(first),'_',timestr,'.xls');
      cd(save_directory);
      xlswrite(m,header,colnames,xlsfilename);
      cd(image_directory);
      end; % i (time) loop
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Exclude the synchoronous movement of the Mito;
    %   
    if (gauss_s(cell_no) > 0),

     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Correction for data collected at 1 hz. Do not correct in data collected after October 5, 2005
        if freq(cell_no) == 1,
            Freq_k = number/(images_per_cycle/2);
        else
            Freq_k = number/images_per_cycle;
        end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        


    	for i = 1:number,
            time1(i) = (i-1)*(1/frequency_hz)/number;
            %time1(i) = (i)*0.32;  
            sinus(i) = sin(2*pi*time1(i)*frequency_hz*Freq_k);
            cosinus(i) = cos(2*pi*time1(i)*frequency_hz*Freq_k);  
            
		end;  
        
        deltaX_backup = deltaX;
        deltaY_backup = deltaY;
        
        for mito_no = 1:length(x_mito),
            
            sumsin  = sum(deltaX_backup(1:number,mito_no).*sinus(:));
            sumcos  = sum(deltaX_backup(1:number,mito_no).*cosinus(:));
            model_x = 2/number*(sumcos*cosinus(:)+sumsin*sinus(:));
            deltaX(1:number,mito_no) = deltaX(1:number,mito_no)- model_x(:);
            
            sumsin  = sum(deltaY_backup(1:number,mito_no).*sinus(:));
            sumcos  = sum(deltaY_backup(1:number,mito_no).*cosinus(:));
            model_y = 2/number*(sumcos*cosinus(:)+sumsin*sinus(:));
            deltaY(1:number,mito_no) = deltaY(1:number,mito_no)- model_y(:);           

        end;
     
%        deltaX
%        deltaX_backup
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %save the spatial positions and corresponding displacements without sine summed over prescribed 'no_of_points'
    %Save the spatial positions and corresponding displacements for ONLY sine, summed over prescribed 'no_of_points'
    %X Y DELTX DELTY
		
%   Obtain deltaX and deltaY summed over prescribed 'no_of_points' to reduce time pts where traction is evaluated

    DX1 = [];
    DY1 = [];        
    DX2 = [];
    DY2 = [];
     m1 = [];
     m2 = [];
     
     		disp('store data as txt-file');
   
 	for i = 1: floor(length(time_buffer)/no_of_points),
                     DX1(i,:)            = sum(deltaX(1:(no_of_points*i),:),1);
                     DY1(i,:)            = sum(deltaY(1:(no_of_points*i),:),1);
                     DX2(i,:)            = sum(deltaX_backup(1:(no_of_points*i),:),1)-DX1(i,:);
                     DY2(i,:)            = sum(deltaY_backup(1:(no_of_points*i),:),1)-DY1(i,:);
                     DX2(i,:)            = DX(i,:)-DX1(i,:);
                     DY2(i,:)            = DY(i,:)-DY1(i,:);
        
        for j = 1:length(x_mito),    
             DispXtxt        = num2str(deltaX(i,j));
             DispYtxt        = num2str(deltaY(i,j));
              m1(j,:)        = [x_mito_selected(j),y_mito_selected(j),(deltaX(i,j)),(deltaY(i,j))];
         end;
  
          header = '';
        colnames = '';

%     Add a leading zero to time points with 4 digits (so that windows arranges files correcly by name)
        if floor(time_buffer(i*no_of_points))/10000 < 1
            timestr = strcat('0',num2str(time_buffer(i*no_of_points)));
        else
            timestr = num2str(time_buffer(i*no_of_points));
       end;

      xlsfilename = strcat('trypsindisp_NoSine_',num2str(first),'_',timestr,'.xls');
      cd(save_directory);
      xlswrite(m1,header,colnames,xlsfilename);
      cd(image_directory);
      end; % i (time) loop
        
        
    end; 
     
  
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %save the Mean Square Displacement (MSD) with time course, Fields:
    %Time   MSD     
    %

    
%     time_span  = [];
%     MSD        = [];
%  
%     for step_num = 1: length(time_buffer),
%         disp(strcat('MSD is being calculated:', num2str(step_num)));
%         time_series = [];
%         MSD_series  = zeros(length(time_buffer)-step_num+1, length(x_mito));
%         for j = 1: (length(time_buffer)-step_num+1),
%             for i = 1:length(x_mito),
%                 sum_deltaX = 0;
%                 sum_deltaY = 0;
%                 for k = j:(j+step_num-1),
%                     sum_deltaX = sum_deltaX + deltaX(k,i);
%                     sum_deltaY = sum_deltaY + deltaY(k,i);
%                 end;            
%                 MSD_series(j,i)  = sum_deltaX.^2 + sum_deltaY.^2;
%             end;
%             time_series = [time_series; (time_buffer(j+step_num-1)-time_buffer(j)+time_buffer(2)-time_buffer(1))*0.1];
%         end;
%         time_span = [time_span;   mean(time_series)];
%         MSD       = [MSD;   mean(MSD_series,1)];
%         
%     end; %for step_num
%         
%         
%     disp('store data as txt-file');
% 	data_file = strcat(save_directory,'TimeCourseMSD_',char(base(cell_no)),'_',num2str(time(m_min)),'_',num2str(time(m_max)),'.csv');
%     fid2 = fopen(data_file,'w');    
%     fprintf(fid2,'%s,%s\n','Time(sec.)','MSD(nm2)'); 
% 
%  	for i = 1: length(time_span),
%         current_line = [];
%         Timetxt = num2str(time_span(i)); %sec.
%         current_line = strcat(current_line,Timetxt,char(44));
%         for j = 1:length(x_mito),    
%             MSDtxt = num2str(MSD(i,j)*1e+6);       %converted to nm;
%             current_line = strcat(current_line,MSDtxt,char(44));
%         end;
%         fprintf(fid2,'%s\n',current_line); %write the new data line;           
% 	end;
%     
%     fclose(fid2);
%     
  
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %save the Mean Square Displacement (MSD) with time course, using Dr. Jim Butler's algorithm. Fields:
    %Time   MSD     
    % 
    
%     time_span  = [];
%     MSD        = [];
%     
%     for step_num = 1: length(time_buffer),
%         disp(strcat('MSD is being calculated:', num2str(step_num)));
%         time_series = [];
%         MSD_series  = [];       
%         for j = 1:step_num:(length(time_buffer)-step_num+1),
%             MSD_row  = zeros(1, length(x_mito));
%             for i = 1:length(x_mito),
%                 sum_deltaX = 0;
%                 sum_deltaY = 0;
%                 for k = j:(j+step_num-1),
%                     sum_deltaX = sum_deltaX + deltaX(k,i);
%                     sum_deltaY = sum_deltaY + deltaY(k,i);
%                 end;            
%                 MSD_row(1,i)  = sum_deltaX.^2 + sum_deltaY.^2;
%             end;
%             time_series = [time_series; (time_buffer(j+step_num-1)-time_buffer(j)+time_buffer(2)-time_buffer(1))*0.1];
%             MSD_series  = [MSD_series; MSD_row];
%         end;
%         time_span = [time_span;   mean(time_series)];
%         MSD       = [MSD;   mean(MSD_series,1)];
%         
%     end; %for step_num
%         
%         
%     disp('store data as txt-file');
% 	data_file = strcat(save_directory,'TimeCourseMSD_NotOverlap_',char(base(cell_no)),'_',num2str(time(m_min)),'_',num2str(time(m_max)),'.csv');
%     fid3 = fopen(data_file,'w');    
%     fprintf(fid3,'%s,%s\n','Time(sec.)','MSD(nm2)'); 
% 
%  	for i = 1: length(time_span),
%         current_line = [];
%         Timetxt = num2str(time_span(i)); %sec.
%         current_line = strcat(current_line,Timetxt,char(44));
%         for j = 1:length(x_mito),    
%             MSDtxt = num2str(MSD(i,j)*1e+6);       %converted to nm;
%             current_line = strcat(current_line,MSDtxt,char(44));
%         end;
%         fprintf(fid3,'%s\n',current_line); %write the new data line;           
% 	end;
%     
%     fclose(fid3);
    
    
    

end; %(cell_no)

end; %(fil)
