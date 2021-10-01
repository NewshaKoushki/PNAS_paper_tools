% Rearrange rawdata to be accessible to cropping
% Check if files are .tif or .bmp
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

fdate = 'rawdata';
folds  = {'Cell0001'}; % cell folder to analyze

main_fold = pwd;
foldAdr = 'D:\Dropbox\Pooya-Rose\Science\TractionForceMicroscopy';

i1 = 700;
i2 = 1499;
j1 = 700;
j2 = 1499;

nchannels = 3;
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

folder1 = 'Mark_and_Find_001';
folder2 = 'Mark_and_Find_002';

for fi1 = 1:length(folds)
    cellfolder	= folds{fi1};
    
    cropinput_Haruka_new;
    
%     d2 = dir(util.gen_addr('..\rawdata'));
    
    cd(util.gen_addr([foldAdr '/' fdate '/' cellfolder '/' folder1]));
    d2a = dir;
    %%
    %============================================================================================
    % Dedrift and then crop
    
    mkdir(strcat([foldAdr '/croppeddata/',cellfolder,'/BF']));
    mkdir(strcat([foldAdr '/croppeddata/',cellfolder,'/BigBF']));
    %%mkdir(strcat('\HTSbart\croppeddata\',cellfolder,'\BigFB'));
    cd([foldAdr '/' fdate '/' cellfolder, '/',folder2]);
    clear d2b d2i d2;
    d2 = dir; si = 1;
    d2b(1)=d2(1);			%d2b is part of d2 and will include [.][..] and original cell files
    d2b(2)=d2(2);
    
    for d2i = 3 : length(d2),
        if ~isempty(findstr(d2(d2i).name,'Position'))
            d2b(si+2)	=d2(d2i);
            d2b(si+2).name 	= d2(d2i).name;
            si = si + 1;
        end;
    end;
    
    
    pos = (si-1)/nchannels;
    
    % Evaluate position by position
    % Change -2 to -3 in line below if there is a .DSStore file
    for posn = 1 : pos,% CHANGED THIS TO GET VALUES FROM FIRST ROW OF DISHES
        % Identifying the trypsin images (Called First)
        
        if  posn <= 9,
            FB = strcat(d2b(pos+2).name(1:27),num2str(posn),'_t9999_ch00.tif');
            BigFB = strcat(d2b(pos+2).name(1:27),num2str(posn),'_t9999_ch01.tif');
            BF = strcat(d2b(pos+2).name(1:27),num2str(posn),'_t9999_ch02.tif');
        else
            FB = strcat(d2b(pos+2).name(1:26),num2str(posn),'_t9999_ch00.tif');
            BigFB = strcat(d2b(pos+2).name(1:26),num2str(posn),'_t9999_ch01.tif');
            BF = strcat(d2b(pos+2).name(1:26),num2str(posn),'_t9999_ch02.tif');
        end;
        
        cd([foldAdr '/' fdate '/' cellfolder '/' folder2]);
        % Dedrift based on top beads
        im1FB 			= double(imread(FB));
        im1FB_S001    	= im1FB(i1:i2,j1:j2);
        
        im1BF 			= double(imread(BF));
        im1BF_S001		= im1BF(i1:i2,j1:j2);
        
        im1BigFB 	    = double(imread(BigFB));
        im1BigFB_S001   = im1BigFB(i1:i2,j1:j2);
        
        
        
        % Identifying the experimental images (Called Second)
        % Change "-2" to "-3" in the line below if there is a .DS Store
        % file
        no_tpts = (length(d2a)-2)/(si-1);
        
        for time = 1 : no_tpts,
            %%            try
            if  posn <= 9,
                if time <= 10,
                    name2 = strcat(d2b(pos+2).name(1:27),num2str(posn),'_t',num2str(time-1,'%04d'),'_ch00.tif');
                else
                    name2 = strcat(d2b(pos+2).name(1:27),num2str(posn),'_t', num2str(time-1,'%04d'),'_ch00.tif');
                    
                end;
            else
                if time <= 10,
                    name2 = strcat(d2b(pos+2).name(1:26),num2str(posn),'_t',num2str(time-1,'%04d'),'_ch00.tif');
                else
                    name2 = strcat(d2b(pos+2).name(1:26),num2str(posn),'_t',num2str(time-1,'%04d'),'_ch00.tif');
                end;
                
            end;
            
            
            cd([foldAdr '/' fdate '/' cellfolder '/' folder1]);
            im2FB 			= double(imread(name2));
            im2a_S001		= im2FB(i1:i2,j1:j2);
            
            cd([foldAdr '/Matlab_programs']);
            % CROP the FB image
            
            % ---------------------------------
            [shx_S001,shy_S001] = disp_on_blocks(im1FB_S001,im2a_S001,size(im1FB_S001,1),0);
            im2b_S001			= im2FB(max(1,i1+shy_S001) : min(size(im2FB,1),i2+shy_S001), max(1,j1+shx_S001) : min(size(im2FB,2),j2+shx_S001));
            % ---------------------------------
            
            if size(im2b_S001,1) == size(im1FB_S001,1) || size(im2b_S001,2) == size(im1FB_S001,2),
                im2FB = im2b_S001;
            else 	im2FB = im2a_S001; disp(' Problem in size! ');
            end;
            
            imwrite(uint8(im1FB_S001),[foldAdr '/croppeddata/' cellfolder '/' FB],'tif');
            imwrite(uint8(im2FB),[foldAdr '/croppeddata/' cellfolder '/' name2],'tif');
            
            % Skip to the next file if displacement has already been computed
            % for the file in question
            
            cd([foldAdr '/displacement/' cellfolder]);
            d2check = dir;
            check = 0;
            for dispcheck = 1 : length(d2check),
                if ~isempty(findstr(d2check(dispcheck).name,name2)),	% name must  contain Pos
                    check = check+1;
                end;
            end;
            
            cd(main_fold);
            if check == 1,
%                 try
                    
                    im1 			= im1FB_S001;
                    im2 			= im2FB;
                    name2          = name2;
                    beads_imcorr_bart;
%                 catch
%                     continue
%                 end;
                
            end;
            
            
            
            cd([foldAdr '/' fdate '/' cellfolder '/' folder1]);
            
            % CROP THE BF IMAGES:
            if  posn <= 9,
                if time <= 10,
                    name4 = strcat(d2b(pos+2).name(1:27),num2str(posn),'_t',num2str(time-1,'%04d'),'_Ch01.tif');
                else
                    name4 = strcat(d2b(pos+2).name(1:27),num2str(posn),'_t',num2str(time-1,'%04d'),'_ch01.tif');
                end;
            else
                if time <= 10,
                    name4 = strcat(d2b(pos+2).name(1:26),num2str(posn),'_t',num2str(time-1,'%04d'),'_Ch01.tif');
                else
                    name4 = strcat(d2b(pos+2).name(1:26),num2str(posn),'_t',num2str(time-1,'%04d'),'_ch01.tif');
                end;
                
            end;
            
            
            im4 			= double(imread(name4));
            im4_S001			= im4(max(1,i1+shy_S001) : min(size(im4,1),i2+shy_S001), ...
                max(1,j1+shx_S001) : min(size(im4,2),j2+shx_S001));
            
            imwrite(uint8(im1BF_S001),[foldAdr '/croppeddata/' cellfolder '/BF/' BF],'tif');
            imwrite(uint8(im4_S001),[foldAdr '/croppeddata/' cellfolder '/BF/' name4],'tif');
            
            
            
            % CROP THE CELL FLUORESCENCE IMAGES:
            if  posn <= 9,
                if time <= 10,
                    name5 = strcat(d2b(pos+2).name(1:27),num2str(posn),'_t',num2str(time-1,'%04d'),'_Ch02.tif');
                else
                    name5 = strcat(d2b(pos+2).name(1:27),num2str(posn),'_t',num2str(time-1,'%04d'),'_Ch02.tif');
                end;
            else
                if time <= 10,
                    name5 = strcat(d2b(pos+2).name(1:26),num2str(posn),'_t',num2str(time-1,'%04d'),'_Ch02.tif');
                else
                    name5 = strcat(d2b(pos+2).name(1:26),num2str(posn),'_t',num2str(time-1,'%04d'),'_Ch02.tif');
                end;
                
            end;
            
            im5 			= double(imread(name5));
            im5_S001	    = im5(max(1,i1+shy_S001) : min(size(im5,1),i2+shy_S001), ...
                max(1,j1+shx_S001) : min(size(im5,2),j2+shx_S001));
            
            
            imwrite(uint8(im1BigFB_S001),[foldAdr '/croppeddata/' cellfolder '/BigBF/' BigFB],'tif');
            imwrite(uint8(im5_S001),[foldAdr '/croppeddata/' cellfolder '/BigBF/' name5],'tif');
            
            
            %%   catch
            %%       continue
            %%       end;
            save(horzcat('workspace','.mat') );
        end;
        
    end;
    
end;
cd(main_fold);
