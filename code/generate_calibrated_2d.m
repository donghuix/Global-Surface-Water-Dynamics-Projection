clear;close all;clc;
[e3sm_input, exportfig] = SetupEnvironment();

domain_file = [e3sm_input 'share/domains/domain.lnd.r05_oEC60to30v3.190418.nc'];
load('../data/par_cal.mat');
xc = ncread(domain_file,'xc');
load('index_lnd.mat');

[m,n] = size(xc);

fc2d = ones(m,n).*0.4;
fover2d = ones(m,n).*0.5;

fc2d(index_lnd)    = fc_cal;
fover2d(index_lnd) = fover_cal;

ncaddvar('../inputdata/surfdata_0.5x0.5_simyr2010_c200624.nc','../inputdata/tmp.nc','pc',fc2d,'NC_DOUBLE');
ncaddvar('../inputdata/tmp.nc','../inputdata/surfdata_0.5x0.5_simyr2010_c200624_surfacewater.nc','fover',fover2d,'NC_DOUBLE');