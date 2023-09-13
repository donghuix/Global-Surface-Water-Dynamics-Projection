clear;close all;clc;

[e3sm_input, exportfig] = SetupEnvironment();

case_names = {'ELMMOS_GLOBE_Surface_Water_847df06a6_cal01.2023-09-08-150838'};

ind = 1;

yr_start = 1974;

elmfiles = dir(['../outputs/' case_names{ind} '/run/*.elm.h0.*nc']);
mosfiles = dir(['../outputs/' case_names{ind} '/run/*.mosart.h0.*nc']);

elmfiles = elmfiles((1999-1974)*12+1:end);
mosfiles = mosfiles((1999-1974)*12+1:end);

for i = 1 : length(elmfiles)
    felm = fullfile(elmfiles(i).folder,elmfiles(i).name);
    fmos = fullfile(mosfiles(i).folder,mosfiles(i).name);
    
    disp([num2str(i) '/' num2str(length(elmfiles))]);
    
    if i == 1
        fh2osfc = ncread(felm,'FH2OSFC');
        flooded = ncread(fmos,'FLOODPLAIN_FRACTION');
    else
        fh2osfc(:,i) = ncread(felm,'FH2OSFC');
        flooded(:,i) = ncread(fmos,'FLOODPLAIN_FRACTION');
    end
end

save(['cal_' num2str(ind) '.mat'],'fh2osfc','flooded');
