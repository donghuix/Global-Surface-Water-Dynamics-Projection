clear;close all;clc;

[e3sm_input, exportfig] = SetupEnvironment();

%case_names = {'ELMMOS_GLOBE_Surface_Water_847df06a6_cal01.2023-09-08-150838'};
case_names = {'GLOBE_Surface_Water_Projection_gfdl-esm4_historical_0b7cdf4760.2023-11-15-211244', ...
              'GLOBE_Surface_Water_Projection_ipsl-cm6a-lr_historical_0b7cdf4760.2023-11-17-211701'};

for icase = 1 : length(case_names)
    strs     = strsplit(case_names{icase},'_');
    model    = strs{5};
    scenario = strs{6};
    
    if strcmp(scenario, 'historical')
        yr_start = 1951;
        yr_read  = 1971;
        yr_end   = 2014;
    elseif strcmp(scenraio, 'ssp126') || strcmp(scenario, 'ssp585')
        yr_start = 2015;
        yr_read  = 2015;
        yr_end   = 2100;
    else
        error([scenario ' is not available!']);
    end
    
    if ~exist(['../data/projection_' model '_' scenario '.mat'],'file')
        elmfiles = dir(['../outputs/' case_names{icase} '/run/*.elm.h0.*nc']);
        mosfiles = dir(['../outputs/' case_names{icase} '/run/*.mosart.h0.*nc']);

        elmfiles = elmfiles((yr_read-yr_start)*12+1:end);
        mosfiles = mosfiles((yr_read-yr_start)*12+1:end);
        assert(length(elmfiles) == length(mosfiles));
        assert(length(elmfiles) == (yr_end - yr_read + 1) * 12);
        
        for i = 1 : length(elmfiles)
            felm = fullfile(elmfiles(i).folder,elmfiles(i).name);
            fmos = fullfile(mosfiles(i).folder,mosfiles(i).name);
            
            disp([model ', ' scenario ': ' num2str(i) '/' num2str(length(elmfiles))]);
            
            if i == 1
                fh2osfc = ncread(felm,'FH2OSFC');
                zwt     = ncread(felm,'ZWT');
                tws     = ncread(felm,'TWS');
                flooded = ncread(fmos,'FLOODPLAIN_FRACTION');
            else
                fh2osfc(:,i) = ncread(felm,'FH2OSFC');
                zwt(:,i)     = ncread(felm,'ZWT');
                tws(:,i)     = ncread(felm,'TWS');
                flooded(:,i) = ncread(fmos,'FLOODPLAIN_FRACTION');
            end
        end
        
        save(['../data/projection_' model '_' scenario '.mat'],'fh2osfc','flooded','zwt','tws');
    end
    
end