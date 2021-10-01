% Rearrange rawdata to be accessible to cropping
% Check if files are .tif or .bmp 
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
fdate = 'rawdata';
folds  = {'Cell99';}; % cell folder to analyze
i1 = 16;
i2 = 975;
j1 = 266;
j2 = 1225;

%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


    for fi1 = 1: length(folds)      
        cellfolder	= folds{fi1}
        
 eval(['cd C:\Singlecell\' fdate '\' cellfolder]); 
  d2 = dir;
  for d2i = 3 : length(d2),
        if d2i-2 <= 9,
       folder = strcat('Mark_and_Find_00',num2str(d2i-2));
        else 
      folder = strcat('Mark_and_Find_0',num2str(d2i-2));      
        end;
       
       eval(['cd C:\Singlecell\' fdate '\' cellfolder '\' folder]); 
       d2a = dir; si = 0; 
       for d2j = 3 : length (d2a)
       %Put a time stamp on every image
         if  ~isempty(findstr(d2a(d2j).name,'_ch')) & ...		% name must contain '_ch'
                ~isempty(findstr(d2a(d2j).name,'.tif')),		%name has to contain '.tif; .bmp for Xinyong'     
 		      name 		= d2a(d2j).name;
              name1     = strcat(name(1:16), '_t',num2str(d2i-2),'.tif'); 
              copyfile (name,strcat('../',name1));
              si = si + 1;
   end;
    end;
    end;
%============================================================================================
% Dedrift and then crop 

    mkdir(strcat('C:\Singlecell\croppeddata\',cellfolder,'\BF'));
%%mkdir(strcat('C:\Singlecell\croppeddata\',cellfolder,'\BigFB'));
    eval(['cd C:\Singlecell\' fdate '\' cellfolder]); 

     
    % Evaluate position by position
    tryp = length(d2)-2;        % The time digit for trypsin images
    for pos = 1 : si/3,% CHANGED THIS TO GET VALUES FROM FIRST ROW OF DISHES
        % Identifying the trypsin images (Called First)

        if  pos <= 9, 
        FB = strcat('Pos00',num2str(pos),'_S001_ch00_t',num2str(tryp),'.tif');
        BF = strcat('Pos00',num2str(pos),'_S001_ch02_t',num2str(tryp),'.tif');
        else
        FB = strcat('Pos0',num2str(pos),'_S001_ch00_t',num2str(tryp),'.tif');
        BF = strcat('Pos0',num2str(pos),'_S001_ch00_t',num2str(tryp),'.tif');            
        end;
        
        %eval(['cd C:\Singlecell\' fdate '\' cellfolder '\' folder2]);
        % Dedrift based on top beads
        im1FB 			= double(imread(FB)); 
        im1FB 			= im1FB(i1:i2,j1:j2);  
        im1BF 			= double(imread(BF)); 
        im1BF 			= im1BF(i1:i2,j1:j2);    
        
          % Identifying the experimental images (Called Second)
        for count = 1 : length(d2)-3, 
        
        if  pos <= 9,     
        name2 = strcat('Pos00',num2str(pos),'_S001_ch00_t',num2str(count),'.tif');     
        else
        name2 = strcat('Pos0',num2str(pos),'_S001_ch00_t',num2str(count),'.tif');    
        end;
        im2FB 			= double(imread(name2)); 
		im2a 			= im2FB(i1:i2,j1:j2);
        
               
        % CROP the FB image

        [shx,shy] 	= disp_on_blocks(im1FB,im2a,size(im1FB,1),0); 
%        [shx,shy] 	= disp_on_blocks(im1FB,im2a,128,0); 
%         shx = 0;
%         shy = 0;
        im2b			= im2FB(max(1,i1+shy) : min(size(im2FB,1),i2+shy), ...
      					max(1,j1+shx) : min(size(im2FB,2),j2+shx));
   if size(im2b,1) == size(im1FB,1) | size(im2b,2) == size(im1FB,2),             
      	im2FB = im2b;
   else 	im2FB = im2a; disp(' Problem in size! ');
   end;
      
     imwrite(uint8(im1FB),['C:\Singlecell\croppeddata\' cellfolder '\' FB],'tif');
     imwrite(uint8(im2FB),['C:\Singlecell\croppeddata\' cellfolder '\' name2],'tif');
     
           
       % CROP THE BF IMAGES:
       if pos <= 9,
           name4 = strcat('Pos00',num2str(pos),'_S001_ch02_t',num2str(count),'.tif');          
       else
       name4 = strcat('Pos0',num2str(pos),'_S001_ch02_t',num2str(count),'.tif');          
       end;
      im4 			= double(imread(name4)); 
   	  im4			= im4(max(1,i1+shy) : min(size(im4,1),i2+shy), ...
      					max(1,j1+shx) : min(size(im4,2),j2+shx));
       imwrite(uint8(im1BF),['C:\Singlecell\croppeddata\' cellfolder '\BF\' BF],'tif');
      imwrite(uint8(im4),['C:\Singlecell\croppeddata\' cellfolder '\BF\' name4],'tif');               
                    
              
% %    Compute the displacements
 try
          im1 			= im1FB; 
          im2 			= im2FB;  
          beads_imcorr_bart;    
 catch
     continue
  end
 
%eval(['cd C:\HTSbart\Displacement' '\' cellfolder]); pwd
   
        end;
   
   
    end;
    end;