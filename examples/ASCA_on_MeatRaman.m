clear; clc; close all; 
load Raman_meat_ASCAdataexample.mat;

D = [temp time pig side]; 
Dlb = {'Temp';'Time';'pig';'side'}; 
tab = checkDesign(D,Dlb); shg; 

M = [1 0 0 0; 0 1 0 0; 1 1 0 0]; 

r = ASCAcat(x,D,Dlb,10,M);

