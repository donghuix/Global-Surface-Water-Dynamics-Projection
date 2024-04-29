clear;close all;clc;

re = 6.37122e6;% Earth radius

[e3sm_input, exportfig] = SetupEnvironment();
xv = ncread('../data/domain_lnd_GLOBE_1d.nc','xv');
yv = ncread('../data/domain_lnd_GLOBE_1d.nc','yv');

fsatpar = load(['../data/par_cal_' num2str(12) '.mat']);
swfpar = load(['../data/par_cal_' num2str(13) '.mat']);

figure;
subplot(2,1,1);
patch(xv,yv,fsatpar.fc_cal,'LineStyle','none'); colorbar;
subplot(2,1,2);
patch(xv,yv,swfpar.fc_cal,'LineStyle','none'); colorbar;

figure;
subplot(2,1,1);
patch(xv,yv,fsatpar.fover_cal,'LineStyle','none'); colorbar; 
subplot(2,1,2);
patch(xv,yv,swfpar.fover_cal,'LineStyle','none'); colorbar;