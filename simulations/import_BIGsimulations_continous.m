clc; 
cd '/Users/mortenarendtrasmussen/Dropbox/Backup/MyDocumentsOnC/MATLAB/mytools/ASCA/simulations/'
clear; 
load BIGsimulation.mat
EXPVARall = EXPVAR;
pvall = pv; 
load BIGsimulation2.mat
ic = sum(isnan(EXPVAR),2)==0 ;
EXPVARall(ic,:) = EXPVAR(ic,:); pvall(ic,:) = pv(ic,:);  
load BIGsimulation3.mat
ic = sum(isnan(EXPVAR),2)==0 ;
EXPVARall(ic,:) = EXPVAR(ic,:); pvall(ic,:) = pv(ic,:);  
load BIGsimulation4.mat
ic = sum(isnan(EXPVAR),2)==0 ;
EXPVARall(ic,:) = EXPVAR(ic,:); pvall(ic,:) = pv(ic,:);  
ic = sum(isnan(EXPVARall),2)==0 ;
sum(ic)

save BIGsimulation_combined.mat EXPVARall pvall DOE



