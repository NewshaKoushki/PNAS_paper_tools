% Input file to read i1,i2,i3,i4 and j1,j2,j3,j4 - Coordinates to crop and
% dedrift for patterened monolayers

   
  switch cellfolder

      case{'Cell0102'}
        
        i1 = 16;
        i2 = 1871;
        j1 = 16;
        j2 = 1871;
        
       case{'Cell0823'}
        
        i1 = 16;
        i2 = 1871;
        j1 = 16;
        j2 = 1871;
 pixelsize     = 0.189; %0.92 for Rama's Inverted 10x ; 1.297 for Xinyong 10x; 0.462 for Bart - 20x 
  young        = 9000; %Pa
  pois         = 0.48;
        
   case{'Cell0824'}
        

        i1=16;
       i2 = 1871;
        j1 = 16;
       j2 = 1871;
 pixelsize     = 0.189; %0.92 for Rama's Inverted 10x ; 1.297 for Xinyong 10x; 0.462 for Bart - 20x 
  young        = 9000; %Pa
  pois         = 0.48;
  
     end;% switch

     

  