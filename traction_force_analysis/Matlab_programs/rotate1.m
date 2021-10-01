function [a_sym,lambda,th] = rotate1(a);

a_sym  = 0.5 * (a + a');

th = atan2((a(1,2) + a(2,1)) , (a(2,2) - a(1,1))) / 2;

U  = [cos(th) sin(th); 
     -sin(th) cos(th)];
   
lambda = U' * a_sym * U;
