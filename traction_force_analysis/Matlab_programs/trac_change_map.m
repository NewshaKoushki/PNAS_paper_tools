%******************************************************************** 
% Program to calculate traction map based on ratio of tractions
% at different time points
% 
% Inputs:
%    Traction maps at different time points
%   
% Outputs:
%    1. % traction map ;
%
% Written by Ramaswamy Krishnan on 08/29/2006
%********************************************************************

clear all; close all;

load tracmap_1501;     %  read data into the my_xy matrix
x = tracmap_1501(:,1);   
