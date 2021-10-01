   %******************************************************************************
   %Changes suggested in fttc_v6.m
   %******************************************************************************

   
   
   %Add this 2 lines at the beginning of the file just after the code line : "young = young * 1e-12;"
   %Please feel free to change this two parameters according to your analysis results;
   MaxAmplitude = 5.0;  %um;
   MaxTraction  = 800;; %Pa;

   
   
   
   
   %Add the following lines around line 170 just after the code: "m1 = max2(sqrt(ux.^2+uy.^2));" 
   set(cbh,'YTick',0:0.5:m1*(m1/MaxAmplitude));
   set(cbh,'YTickLabel',{str2num(num2str(0:0.5*(MaxAmplitude/m1):m1,2))});  
   CLim      = get(gca,'CLim');    
   set(gca,'CLim',[CLim(1) CLim(2)*MaxAmplitude/m1]);
   
   
   
   
   
   
   
   
   
   %Add the following lines around line 185 just after the code: "m2 = max2(sqrt(tx.^2+ty.^2)*1e12);"
   set(cbh,'YTick',0:40:m2*(m2/MaxTraction));
   set(cbh,'YTickLabel',{round(0:40*(MaxTraction/m2):m2)});  
   CLim      = get(gca,'CLim');    
   set(gca,'CLim',[CLim(1) CLim(2)*MaxTraction/m2]);
  
  
  
  
  
  
  
   %Add the following lines around line 549 just after the code: "m3 = max2(sqrt(tx.^2+ty.^2)*1e12);"
   set(cbh,'YTick',0:40:m3*(m3/MaxTraction));
   set(cbh,'YTickLabel',{round(0:40*(MaxTraction/m3):m3)});  
   CLim      = get(gca,'CLim');    
   set(gca,'CLim',[CLim(1) CLim(2)*MaxTraction/m3]);
   
         