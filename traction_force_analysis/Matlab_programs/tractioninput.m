
pois = 0.48; % do not change
subcellnum = 1; % do not change

    
  switch cellfolder

   case{'Cell1'}
% 10/14/09 Cell1 of Daphne's experiments (Date:; Condition: ) 
     first          = 'trypsin';   								% Reference image (tryp) don't change
     i1 			= 200; 	% Top  add 1 to pixel point from paint
     i2 			= 903; 	% Bottom
     j1 			= 300; 	% Left  add 1 to pixel point from paint
     j2 			= 1003;	% Right
   cell_img         = {'crop_phase.tif';}; % don't change
   RMSTmin          = 10; % don't change
    
 
     end;% switch
