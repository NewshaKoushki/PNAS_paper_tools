pixelsize = 0.163;
young = 11000;
pois = 0.48;
young = young * 1e-12;
subcellnum = 1;

    
  switch cellfolder
      
case{'Cell12'}

    i1 		= 1; 	% Top  add 1 to pixel point from paint
    i2 			= 600; 	% Bottom
    j1 			= 1; 	% Left  add 1 to pixel point from paint
    j2 			= 800;	% Right
 min_time = [ 08558;  09158;  09598;        ];
 max_time = [ 08629;  09260;  09700;            ];
step     = [ 1;   1;  1;            ];
gauss_s  = [ 50;  50;   50              ];
freq_s   = [0.75;   0.75;   0.75  ];
exposure_time = [0.1;0.1; 0.1         ];   
event    = [1;  1; 1 ];


case{'Cell2'}
    % Amplitude = 75 Gauss    
%    first   	= [10771];   								% Reference image (tryp) 
%    second  	= [07984 08736 12113 17623 20897]; % b, h5-1nd, col
    i1 			= 1; 	% Top  add 1 to pixel point from paint
    i2 			= 600; 	% Bottom
    j1 			= 1; 	% Left  add 1 to pixel point from paint
    j2 			= 600;	% Right
   min_time     = [ 01629; 02029; 03375; 05164];
   max_time     = [ 01667; 02069; 03415; 05262];    
   gauss_s  = [ 75; 75;75;  75  ];
   freq_s   = [0.75; 0.75; 0.75; 0.75  ];
   step     = [1; 1; 1; 1];
exposure_time = [0.1; 0.1; 0.1; 0.1];   
event    = [1; 1;1; 1];

case{'Cell3'}
    % Amplitude = 75 Gauss    
%    first   	= [10771];   								% Reference image (tryp) 
%    second  	= [07984 08736 12113 17623 20897]; % b, h5-1nd, col
    i1 			= 1; 	% Top  add 1 to pixel point from paint
    i2 			= 600; 	% Bottom
    j1 			= 1; 	% Left  add 1 to pixel point from paint
    j2 			= 600;	% Right
   min_time     = [ 01929; 02276; 03076; 03436; 05391];
   max_time     = [ 01969; 02316; 03116; 03476; 05431];    
   gauss_s  = [ 15; 15;15;  15; 15  ];
   freq_s   = [0.75; 0.75; 0.75; 0.75; 0.75  ];
   step     = [1; 1; 1; 1; 1];
exposure_time = [0.1; 0.1; 0.1; 0.1; 0.1];   
event    = [1; 1;1; 1; 1];


case{'Cell4'}
    % Amplitude = 75 Gauss    
%    first   	= [10771];   								% Reference image (tryp) 
%    second  	= [07984 08736 12113 17623 20897]; % b, h5-1nd, col
    i1 			= 1; 	% Top  add 1 to pixel point from paint
    i2 			= 600; 	% Bottom
    j1 			= 1; 	% Left  add 1 to pixel point from paint
    j2 			= 600;	% Right
   min_time     = [ 01656; 02030; 05180];
   max_time     = [ 01710; 02083; 05247];    
   gauss_s  = [ 75; 75;75;    ];
   freq_s   = [0.75; 0.75; 0.75  ];
   step     = [1; 1; 1];
exposure_time = [0.1; 0.1; 0.1];   
event    = [1; 1;1];
   
case{'Cell5'}
        i1 		= 1; 	% Top  add 1 to pixel point from paint
    i2 			= 600; 	% Bottom
    j1 			= 1; 	% Left  add 1 to pixel point from paint
    j2 			= 600;	% Right
 min_time = [ 08464;    09036;  09662;  12978            ];
 max_time = [ 08534;    09131;  09743;  13046            ];
step     = [ 1;   1;  1;1              ];
gauss_s  = [ 75;  75;    75;75              ];
freq_s   = [0.75;   0.75;   0.75;0.75  ];
exposure_time = [0.1;0.1; 0.1;0.1         ];   
event    = [1;  1; 1;1  ];

case{'Cell7'}
      i1 		= 1; 	% Top  add 1 to pixel point from paint
    i2 			= 600; 	% Bottom
    j1 			= 1; 	% Left  add 1 to pixel point from paint
    j2 			= 600;	% Right
 min_time = [ 04705;    05305;  05785;  07974            ];
 max_time = [ 04781;    05387;  05867;  08055            ];
step     = [ 1;   1;  1;1              ];
gauss_s  = [ 15;  15;    15;15              ];
freq_s   = [0.75;   0.75;   0.75;0.75  ];
exposure_time = [0.1;0.1; 0.1;0.1         ];   
event    = [1;  1; 1;1  ];

case{'Cell9'}
      i1 		= 1; 	% Top  add 1 to pixel point from paint
    i2 			= 600; 	% Bottom
    j1 			= 1; 	% Left  add 1 to pixel point from paint
    j2 			= 600;	% Right
 min_time = [ 01808;    02167;  02821; 03421; 05238  ];
 max_time = [ 01846;    02207;  02861; 03461; 05278  ];
step     = [ 1;   1;  1; 1; 1;    ];
gauss_s  = [ 15;  15;   15     ;  15; 15     ];
freq_s   = [0.75;   0.75;0.75; 0.75; 0.75];
exposure_time = [0.1;0.1;0.1 ; 0.1; 0.1     ];   
event    = [1;  1; 1 ; 1; 1];
    
case{'Cell12'}
     i1 		= 1; 	% Top  add 1 to pixel point from paint
    i2 			= 600; 	% Bottom
    j1 			= 1; 	% Left  add 1 to pixel point from paint
    j2 			= 600;	% Right
 min_time = [ 09073;    09676;  10237; 12485;      ];
 max_time = [ 09142;    09745;  10345; 12570;      ];
step     = [ 1;   1;  1;   1   ];
gauss_s  = [ 75;  75;   75; 75             ];
freq_s   = [0.75;   0.75;0.75; 0.75  ];
exposure_time = [0.1;0.1;0.1; 0.1         ];   
event    = [1;  1; 1; 1 ];
     
    
case{'Cell13'}
     i1 		= 1; 	% Top  add 1 to pixel point from paint
    i2 			= 600; 	% Bottom
    j1 			= 1; 	% Left  add 1 to pixel point from paint
    j2 			= 600;	% Right
 min_time = [ 03252; 03746; 04320;                          06736;                        ];
 max_time = [ 03320; 03827; 04401;                          06831;                        ];
step     = [ 1;         1; 1            ;         1;                             ];
gauss_s  = [ 75;                              75;      75;75                       ];
freq_s   = [0.75;                              0.75;     0.75;0.75                      ];
exposure_time = [0.1;                       0.1;             0.1;0.1             ];   
event    = [1;  1;  1; 1 ];

 case{'Cell20'}
% Amplitude = 75 Gauss    
    first   	= [13575];   								% Reference image (tryp) 
%    second  	= [07984 08736 12113 17623 20897]; % b, h5-1nd, col
    i1 			= 1; 	% Top  add 1 to pixel point from paint
    i2 			= 600; 	% Bottom
    j1 			= 1; 	% Left  add 1 to pixel point from paint
    j2 			= 600;	% Right
   min_time     = [ 05300; 05673; 05727; 06927; ];
   max_time     = [ 05328; 05700; 05753; 06953; ];    
   gauss_s  = [ 75;  75; 75; 75   ];
   freq_s   = [0.75;  0.75; 0.75; 0.75 ];
   step     = [1; 1; 1; 1];
exposure_time = [0.1; 0.1; 0.1; 0.1];   
event    = [1; 1; 1; 1];
   
 case{'Cell21'}
% Amplitude = 75 Gauss    
    first   	= [14976];   								% Reference image (tryp) 
%    second  	= [07984 08736 12113 17623 20897]; % b, h5-1nd, col
    i1 			= 1; 	% Top  add 1 to pixel point from paint
    i2 			= 600; 	% Bottom
    j1 			= 1; 	% Left  add 1 to pixel point from paint
    j2 			= 600;	% Right
   min_time     = [ 04293; ];
   max_time     = [ 04320; ];    
   gauss_s  = [ 75;     ];
   freq_s   = [0.75;   ];
   step     = [1; ];
exposure_time = [0.1; ];   
event    = [1; ];

case{'Cell26'}
    % Amplitude = 75 Gauss    
%    first   	= [10771];   								% Reference image (tryp) 
%    second  	= [07984 08736 12113 17623 20897]; % b, h5-1nd, col
    i1 			= 1; 	% Top  add 1 to pixel point from paint
    i2 			= 600; 	% Bottom
    j1 			= 1; 	% Left  add 1 to pixel point from paint
    j2 			= 600;	% Right
   min_time     = [ 02197; 02542; 03063; 03690; 05402];
   max_time     = [ 02237; 02583; 03103; 03730; 05442];    
   gauss_s  = [ 75; 75;75;  75; 75  ];
   freq_s   = [0.75; 0.75; 0.75; 0.75; 0.75  ];
   step     = [1; 1; 1; 1; 1];
exposure_time = [0.1; 0.1; 0.1; 0.1; 0.1];   
event    = [1; 1;1; 1; 1];

case{'Cell35'}
    % Amplitude = 75 Gauss    
%    first   	= [10771];   								% Reference image (tryp) 
%    second  	= [07984 08736 12113 17623 20897]; % b, h5-1nd, col
    i1 			= 1; 	% Top  add 1 to pixel point from paint
    i2 			= 600; 	% Bottom
    j1 			= 1; 	% Left  add 1 to pixel point from paint
    j2 			= 600;	% Right
   min_time     = [ 01694; 02080; 02827; 03454; 05603];
   max_time     = [ 01734; 02120; 02867; 03494; 05643];    
   gauss_s  = [ 75; 75;75;  75; 75  ];
   freq_s   = [0.75; 0.75; 0.75; 0.75; 0.75  ];
   step     = [1; 1; 1; 1; 1];
exposure_time = [0.1; 0.1; 0.1; 0.1; 0.1];   
event    = [1; 1;1; 1; 1];

end;% switch

