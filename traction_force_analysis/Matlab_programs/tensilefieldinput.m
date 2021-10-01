pixelsize = 0.3263; %0.3263 for 20X zoom; 0.163 for 40X zoom
young = 4000;
pois = 0.48;
young = young * 1e-12;
subcellnum = 1;

    
  switch cellfolder

   case{'Cell1'}
% 5/1/06 Cell1 ; 75 Gaussp
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;

      
  case{'Cell2'}
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
     case{'Cell3'}
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
        case{'Cell4'}
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;

           case{'Cell5'}
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;

              case{'Cell6'}
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
              case{'Cell7'}
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
              case{'Cell8'}
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
                 case{'Cell9'}
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
                    case{'Cell10'}
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
                       case{'Cell11'}
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
                          case{'Cell12'}
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
                           case{'Cell13'}
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
                           case{'Cell14'}
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
                           case{'Cell15'}
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
                           case{'Cell16'}
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
                           case{'Cell17'}
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
                           case{'Cell18'}
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
                           case{'Cell19'}
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
                           case{'Cell20'}
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
                           case{'Cell21'}
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
                           case{'Cell22'}
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
                              case{'Cell23'}
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
                                 case{'Cell24'}
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
                                 case{'Cell25'}
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
                                 case{'Cell26'}
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
                                 case{'Cell27'}
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
                                 case{'Cell28'}
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
                                 case{'Cell29'}
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
                                 case{'Cell30'}
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
                                 case{'Cell31'}
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
                                 case{'Cell32'}
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
                                 case{'Cell33'}
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
                                 case{'Cell34'}
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
                                 case{'Cell35'}
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
                                 case{'Cell36'}
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
                                 case{'Cell37'}
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
                                 case{'Cell38'}
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
   
                                    case{'Cell39'}
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
                                    case{'Cell40'}
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
                                    case{'Cell41'}
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
   
                                    case{'Cell42'}
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
   
                                    case{'Cell43'}
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
   
                                    case{'Cell44'}
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
                                          case{'Cell45'}
% 5/1/06 Cell1 ; 75 Gaussp
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
                                        case{'Cell46'}
% 5/1/06 Cell1 ; 75 Gaussp
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
                                        case{'Cell47'}
% 5/1/06 Cell1 ; 75 Gaussp
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
                                        case{'Cell48'}
% 5/1/06 Cell1 ; 75 Gaussp
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
                                        case{'Cell49'}
% 5/1/06 Cell1 ; 75 Gaussp
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
                                        case{'Cell50'}
% 5/1/06 Cell1 ; 75 Gaussp
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
                                        case{'Cell51'}
% 5/1/06 Cell1 ; 75 Gaussp
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
                                        case{'Cell52'}
% 5/1/06 Cell1 ; 75 Gaussp
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
                                        case{'Cell53'}
% 5/1/06 Cell1 ; 75 Gaussp
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
                                        case{'Cell54'}
% 5/1/06 Cell1 ; 75 Gaussp
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
                                        case{'Cell55'}
% 5/1/06 Cell1 ; 75 Gaussp
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
                                        case{'Cell56'}
% 5/1/06 Cell1 ; 75 Gaussp
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
                                        case{'Cell57'}
% 5/1/06 Cell1 ; 75 Gaussp
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
                                        case{'Cell58'}
% 5/1/06 Cell1 ; 75 Gaussp
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
                                        case{'Cell59'}
% 5/1/06 Cell1 ; 75 Gaussp
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
                                        case{'Cell60'}
% 5/1/06 Cell1 ; 75 Gaussp
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
                                        case{'Cell61'}
% 5/1/06 Cell1 ; 75 Gaussp
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
                                        case{'Cell62'}
% 5/1/06 Cell1 ; 75 Gaussp
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
                                        case{'Cell63'}
% 5/1/06 Cell1 ; 75 Gaussp
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
                                        case{'Cell64'}
% 5/1/06 Cell1 ; 75 Gaussp
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
                                        case{'Cell65'}
% 5/1/06 Cell1 ; 75 Gaussp
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
                                        case{'Cell66'}
% 5/1/06 Cell1 ; 75 Gaussp
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
                                        case{'Cell67'}
% 5/1/06 Cell1 ; 75 Gaussp
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
                                        case{'Cell68'}
% 5/1/06 Cell1 ; 75 Gaussp
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
                                        case{'Cell69'}
% 5/1/06 Cell1 ; 75 Gaussp
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
                                        case{'Cell70'}
% 5/1/06 Cell1 ; 75 Gaussp
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
                                        case{'Cell71'}
% 5/1/06 Cell1 ; 75 Gaussp
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
                                        case{'Cell72'}
% 5/1/06 Cell1 ; 75 Gaussp
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
   
                                        case{'Cell73'}
% 5/1/06 Cell1 ; 75 Gaussp
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 16; 	% Top  add 1 to pixel point from paint
     i2 			= 783; 	% Bottom
     j1 			= 16; 	% Left  add 1 to pixel point from paint
     j2 			= 783;	% Right
   min_time      = [ 00200;];
   max_time      = [ 20000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
  
                                         case{'Cell120'}
% 5/1/06 Cell1 ; 75 Gaussp
     first       = [99999];   								% Reference image (tryp) 
%     second    	= [07908 08672 11802 15205 18406]; % b, h5-1nd, col
     i1 			= 212; 	% Top  add 1 to pixel point from paint
     i2 			= 979; 	% Bottom
     j1 			= 542; 	% Left  add 1 to pixel point from paint
     j2 			= 1309;	% Right
   min_time      = [ 00200;];
   max_time      = [ 40000;];   
   cell_img      = {'crop_phase.tif';};
   freq          = [ 0;];
   RMSTmin       = 10;
  young = 11000 * 1e-12;
end;% switch
