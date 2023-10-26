clear;close all;clc;

[e3sm_input, exportfig] = SetupEnvironment();

%case_names = {'ELMMOS_GLOBE_Surface_Water_847df06a6_cal01.2023-09-08-150838'};
case_names = {'ELMMOS_GLOBE_Surface_Water_847df06_cal01_compy.2023-09-13-144635',...
              'ELMMOS_GLOBE_Surface_Water_847df06_cal02_compy.2023-09-17-214204',...
              'ELMMOS_GLOBE_Surface_Water_847df06_cal03_compy.2023-09-18-101542',...
              'ELMMOS_GLOBE_Surface_Water_847df06_cal04_compy.2023-10-02-093842',...
              'ELMMOS_GLOBE_Surface_Water_847df06_cal05_compy.2023-10-09-110846',...
              'ELMMOS_GLOBE_Surface_Water_847df06_cal06_compy.2023-10-10-082930',...
              'ELMMOS_GLOBE_Surface_Water_847df06_cal07_compy.2023-10-15-095218',...
              'ELMMOS_GLOBE_Surface_Water_847df06_cal08_compy.2023-10-15-122907'};
yr_start = 1973;
yr_read  = 1993;

for icase = 1 : length(case_names)
    
    if ~exist(['../data/cal_' num2str(icase) '.mat'],'file')
        elmfiles = dir(['../outputs/' case_names{icase} '/run/*.elm.h0.*nc']);
        mosfiles = dir(['../outputs/' case_names{icase} '/run/*.mosart.h0.*nc']);
        
        elmfiles = elmfiles((yr_read-yr_start)*12+1:end);
        mosfiles = mosfiles((yr_read-yr_start)*12+1:end);
        
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
        
        save(['../data/cal_' num2str(icase) '.mat'],'fh2osfc','flooded');
    end
end
