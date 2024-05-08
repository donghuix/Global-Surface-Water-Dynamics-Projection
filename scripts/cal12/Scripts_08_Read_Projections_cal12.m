clear;close all;clc;

addpath('../../code');
[e3sm_input, exportfig] = SetupEnvironment();

case_names = {'GLOBE_Surface_Water_Projection_cal12_gfdl-esm4_historical_0b7cdf4760.2024-04-04-145707',     ...
              'GLOBE_Surface_Water_Projection_cal12_gfdl-esm4_ssp126_0b7cdf4760.2024-04-06-163000',         ...
              'GLOBE_Surface_Water_Projection_cal12_gfdl-esm4_ssp585_0b7cdf4760.2024-04-06-163358'};

for icase = 1 : length(case_names)
    strs     = strsplit(case_names{icase},'_');
    model    = strs{6};
    scenario = strs{7};
    
    if strcmp(scenario, 'historical')
        yr_start = 1951;
        yr_read  = 1971;
        yr_end   = 2014;
    elseif strcmp(scenario, 'ssp126') || strcmp(scenario, 'ssp585')
        yr_start = 2015;
        yr_read  = 2015;
        yr_end   = 2100;
    else
        error([scenario ' is not available!']);
    end
    
    elmfiles = dir(['./outputs/' case_names{icase} '/run/*.elm.h0.*nc']);
    mosfiles = dir(['./outputs/' case_names{icase} '/run/*.mosart.h0.*nc']);

    elmfiles = elmfiles((yr_read-yr_start)*12+1:end);
    mosfiles = mosfiles((yr_read-yr_start)*12+1:end);
    assert(length(elmfiles) == length(mosfiles));
    assert(length(elmfiles) == (yr_end - yr_read + 1) * 12);
        
    if ~exist(['./projection_cal12_' model '_' scenario '.mat'],'file')
        
        for i = 1 : length(elmfiles)
            felm = fullfile(elmfiles(i).folder,elmfiles(i).name);
            fmos = fullfile(mosfiles(i).folder,mosfiles(i).name);
            
            disp([model ', ' scenario ': ' num2str(i) '/' num2str(length(elmfiles))]);
            
            if i == 1
                fh2osfc = ncread(felm,'FH2OSFC');
                zwt     = ncread(felm,'ZWT');
                tws     = ncread(felm,'TWS');
                rain    = ncread(felm,'RAIN');
                tsa     = ncread(felm,'TSA');
                fsat    = ncread(felm,'FSAT');
                perch   = ncread(felm,'ZWT_PERCH');
                flooded = ncread(fmos,'FLOODPLAIN_FRACTION');
            else
                fh2osfc(:,i) = ncread(felm,'FH2OSFC');
                zwt(:,i)     = ncread(felm,'ZWT');
                tws(:,i)     = ncread(felm,'TWS');
                rain(:,i)    = ncread(felm,'RAIN');
                tsa(:,i)     = ncread(felm,'TSA');
                fsat(:,i)    = ncread(felm,'FSAT');
                perch(:,i)   = ncread(felm,'ZWT_PERCH');
                flooded(:,i) = ncread(fmos,'FLOODPLAIN_FRACTION');
            end
        end
        
        save(['./projection_cal12_' model '_' scenario '.mat'],'fh2osfc','flooded','zwt','tws','rain','tsa','fsat','perch');
    end
    
    if ~exist(['./projection_cal12_' model '_' scenario '_ice.mat'],'file')
    
        for i = 1 : length(elmfiles)
            felm = fullfile(elmfiles(i).folder,elmfiles(i).name);
            fmos = fullfile(mosfiles(i).folder,mosfiles(i).name);
            
            disp([model ', ' scenario ': ' num2str(i) '/' num2str(length(elmfiles))]);
            
            if i == 1
                qsnow   = ncread(felm,'QSNOMELT');
                frost   = ncread(felm,'FROST_TABLE');
                fsno    = ncread(felm,'FSNO');
                snow    = ncread(felm,'SNOW');
                tmp     = ncread(felm,'SOILICE');
                soilice = nanmean(tmp,2);
                surfice = tmp(:,1);
            else
                qsnow(:,i)   = ncread(felm,'QSNOMELT');
                frost(:,i)   = ncread(felm,'FROST_TABLE');
                fsno(:,i)    = ncread(felm,'FSNO');
                snow(:,i)    = ncread(felm,'SNOW');
                tmp          = ncread(felm,'SOILICE');
                soilice(:,i) = nanmean(tmp,2);
                surfice(:,i) = tmp(:,1);
            end
        end
        
        save(['./projection_cal12_' model '_' scenario '_ice.mat'],'qsnow','frost','fsno','snow','soilice','surfice');
        
    end
    
    if ~exist(['./projection_cal12_' model '_' scenario '_flux.mat'],'file')
        
        for i = 1 : length(elmfiles)
            felm = fullfile(elmfiles(i).folder,elmfiles(i).name);
            fmos = fullfile(mosfiles(i).folder,mosfiles(i).name);
            
            disp([model ', ' scenario ': ' num2str(i) '/' num2str(length(elmfiles))]);
            
            if i == 1
                latent      = ncread(felm,'EFLX_LH_TOT');
                fgr         = ncread(felm,'FGR');
                tmp         = ncread(felm,'TSOI');
                soiltemp    = nanmean(tmp,2);
                surftsoi    = tmp(:,1);
            else
                latent(:,i) = ncread(felm,'EFLX_LH_TOT');
                fgr(:,i)    = ncread(felm,'FGR');
                tmp         = ncread(felm,'TSOI');
                soiltemp(:,i) = nanmean(tmp,2);
                surftsoi(:,i) = tmp(:,1);
            end
        end
        
        save(['./projection_cal12_' model '_' scenario '_flux.mat'],'latent','fgr','soiltemp','surftsoi');
    
    end
    
    if ~exist(['./projection_cal12_' model '_' scenario '_energy.mat'],'file')
        
        for i = 1 : length(elmfiles)
            felm = fullfile(elmfiles(i).folder,elmfiles(i).name);
            fmos = fullfile(mosfiles(i).folder,mosfiles(i).name);
            
            disp([model ', ' scenario ': ' num2str(i) '/' num2str(length(elmfiles))]);
            
            if i == 1
                fsds      = ncread(felm,'FSDS');
                fsa       = ncread(felm,'FSA');
                fsr       = ncread(felm,'FSR');
                fire      = ncread(felm,'FIRE');
            else
                fsds(:,i) = ncread(felm,'FSDS');
                fsa(:,i)  = ncread(felm,'FSA');
                fsr(:,i)  = ncread(felm,'FSR');
                fire(:,i) = ncread(felm,'FIRE');

            end
        end
        
        save(['./projection_cal12_' model '_' scenario '_energy.mat'],'fsds','fsa','fsr','fire');
    
    end
    
end