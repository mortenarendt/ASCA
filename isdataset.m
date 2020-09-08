function out = isdataset(in)
%ISDATASET Test for dataset, returns true if 'in' is a PLStoolbox dataset
%(Eigenvector Research). 
%
%I/O: out = isdataset(in);

%Copyright Eigenvector Research 2006-2019
%Licensee shall not re-compile, translate or convert "M-files" contained
% in PLS_Toolbox for use with any software other than MATLAB�, without
% written permission from Eigenvector Research, Inc.

out = isa(in,'dataset');
