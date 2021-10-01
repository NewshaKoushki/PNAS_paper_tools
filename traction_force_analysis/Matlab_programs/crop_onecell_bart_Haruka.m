% Rearrange rawdata to be accessible to cropping
% Check if files are .tif or .bmp 
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
fdate = 'rawdata';
folds  = {'Cell0518'}; % cell folder to analyze 

nchannels = 3;
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


    for fi1 = 1: length(folds)      
        cellfolder	= folds{fi1}
        
        
     eval(['cropinput']);
        
 eval(['cd C:\Singlecell\' fdate '\' cellfolder]); 
  d2 = dir;,
       folder1 = strcat('Mark_and_Find_001');
       folder2 = strcat('Mark_and_Find_002');
         
       eval(['cd \Singlecell\' fdate '\' cellfolder '\' folder1]); 
       d2a = dir; si = 0;
    
%============================================================================================
% Dedrift and then crop 

    mkdir(strcat('C:\Singlecell\croppeddata\',cellfolder,'\BF'));
    mkdir(strcat('C:\Singlecell\croppeddata\',cellfolder,'\BigFB'));
    %%mkdir(strcat('\HTSbart\croppeddata\',cellfolder,'\BigFB'));
    eval(['cd C:\Singlecell\' fdate '\' cellfolder, '\',folder2]); 
    d2tryp = dir;
    si = 0,
    for count = 1 : length(d2tryp),
    if  ~isempty(findstr(d2tryp(count).name,'Pos'))&...
            ~isempty(findstr(d2tryp(count).name,'t99')),% name must  contain Pos
          	si = si + 1;
      end;    
    end;
    pos = si/nchannels
     
    % Evaluate position by position
    % Change -2 to -3 in line below if there is a .DSStore file
    for posn = 1 : pos,% CHANGED THIS TO GET VALUES FROM FIRST ROW OF DISHES
        % Identifying the trypsin images (Called First)

        if  posn <= 9, 
        FB = strcat('Pos00',num2str(posn),'_S001_t99_ch00.tif');
        BigFB = strcat('Pos00',num2str(posn),'_S001_t99_ch01.tif');
        BF = strcat('Pos00',num2str(posn),'_S001_t99_ch02.tif');
        else
        FB = strcat('Pos0',num2str(posn),'_S001_t99_ch00.tif');
        BigFB = strcat('Pos0',num2str(posn),'_S001_t99_ch01.tif');
        BF = strcat('Pos0',num2str(posn),'_S001_t99_ch02.tif');
        end;
        
        eval(['cd C:\Singlecell\' fdate '\' cellfolder '\' folder2]);
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
        no_tpts = (length(d2a)-2)/si
        for time = 1 : no_tpts,    
%%            try
        if  posn <= 9,
                  if time <= 10,
        name2 = strcat('Pos00',num2str(posn),'_S001_t0',num2str(time-1),'_ch00.tif')
                  else
                      name2 = strcat('Pos00',num2str(posn),'_S001_t',num2str(time-1),'_ch00.tif')
                  end;
        else
            if time <= 10,
        name2 = strcat('Pos0',num2str(posn),'_S001_t0',num2str(time-1),'_ch00.tif')
                  else
                      name2 = strcat('Pos0',num2str(posn),'_S001_t',num2str(time-1),'_ch00.tif')
                  end;
        
        end;
        
                    
        eval(['cd C:\Singlecell\' fdate '\' cellfolder '\' folder1]);
        im2FB 			= double(imread(name2)); 
		im2a_S001		= im2FB(i1:i2,j1:j2);

        
   
        
        % CROP the FB image

        [shx_S001,shy_S001] 	= disp_on_blocks(im1FB_S001,im2a_S001,size(im1FB_S001,1),0); 
 
%        [shx,shy] 	= disp_on_blocks(im1FB,im2a,128,0); 
%         shx = 0;
%         shy = 0;
        im2b_S001			= im2FB(max(1,i1+shy_S001) : min(size(im2FB,1),i2+shy_S001), ...
      					max(1,j1+shx_S001) : min(size(im2FB,2),j2+shx_S001));
                    

        
   if size(im2b_S001,1) == size(im1FB_S001,1) | size(im2b_S001,2) == size(im1FB_S001,2),             
      	im2FB = im2b_S001;
   else 	im2FB = im2a_S001; disp(' Problem in size! ');
   end;
   
     imwrite(uint8(im1FB_S001),['C:\Singlecell\croppeddata\' cellfolder '\' FB],'tif');
     imwrite(uint8(im2FB),['C:\Singlecell\croppeddata\' cellfolder '\' name2],'tif');
     
% Skip to the next file if displacement has already been computed
% for the file in question
        
      eval(['cd C:\Singlecell\displacement\' cellfolder]);
      d2check = dir; 
      check = 0;
             for dispcheck = 1 : length(d2check),
            if ~isempty(findstr(d2check(dispcheck).name,name2)),	% name must  contain Pos
            check = check+1;
            end;
       end;
       
        
       if check == 1,
         try  
     
         im1 			= im1FB_S001; 
         im2 			= im2FB;
         name2          = name2;
         beads_imcorr_bart;    
         catch
              continue
        end;
       
       end;
         
  
       
       eval(['cd C:\Singlecell\' fdate '\' cellfolder '\' folder1]);
           
       % CROP THE BF IMAGES:
       if posn <= 9,
        name4 = strcat('Pos00',num2str(posn),'_S001_t0',num2str(time-1),'_ch01.tif');
           else
        name4 = strcat('Pos0',num2str(posn),'_S001_t0',num2str(time-1),'_ch01.tif'); 
       end;
       
      im4 			= double(imread(name4)); 
   	  im4_S001			= im4(max(1,i1+shy_S001) : min(size(im4,1),i2+shy_S001), ...
      					max(1,j1+shx_S001) : min(size(im4,2),j2+shx_S001));
     
      imwrite(uint8(im1BF_S001),['C:\Singlecell\croppeddata\' cellfolder '\BF\' BF],'tif');
      imwrite(uint8(im4_S001),['C:\Singlecell\croppeddata\' cellfolder '\BF\' name4],'tif');

      
           
       % CROP THE CELL FLUORESCENCE IMAGES:
       if posn <= 9,
        name5 = strcat('Pos00',num2str(posn),'_S001_t0',num2str(time-1),'_ch02.tif');
           else
        name5 = strcat('Pos0',num2str(posn),'_S001_t0',num2str(time-1),'_ch02.tif'); 
       end;
       
      im5 			= double(imread(name5)); 
   	  im5_S001	    = im5(max(1,i1+shy_S001) : min(size(im5,1),i2+shy_S001), ...
      					max(1,j1+shx_S001) : min(size(im5,2),j2+shx_S001));
     
      
      imwrite(uint8(im1BigFB_S001),['C:\Singlecell\croppeddata\' cellfolder '\BigFB\' BigFB],'tif');
      imwrite(uint8(im5_S001),['C:\Singlecell\croppeddata\' cellfolder '\BigFB\' name5],'tif');
     
   
%%   catch
%%       continue
%%       end;
        end;
        end;
   
    end;