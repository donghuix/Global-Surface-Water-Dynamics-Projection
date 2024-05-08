clear;close all;clc;

[e3sm_input, exportfig] = SetupEnvironment();

%case_names = {'ELMMOS_GLOBE_Surface_Water_847df06a6_cal01.2023-09-08-150838'};
case_names = {'GLOBE_Surface_Water_Projection_gfdl-esm4_historical_0b7cdf4760.2024-02-22-134555',     ...
              'GLOBE_Surface_Water_Projection_gfdl-esm4_ssp126_0b7cdf4760.2024-02-23-081731',         ...
              'GLOBE_Surface_Water_Projection_gfdl-esm4_ssp585_0b7cdf4760.2024-02-23-082732',         ...
              'GLOBE_Surface_Water_Projection_ipsl-cm6a-lr_historical_0b7cdf4760.2023-11-17-211701',  ...
              'GLOBE_Surface_Water_Projection_ipsl-cm6a-lr_ssp126_0b7cdf4760.2024-02-14-110215',      ...
              'GLOBE_Surface_Water_Projection_ipsl-cm6a-lr_ssp585_0b7cdf4760.2024-02-14-110831',      ...
              'GLOBE_Surface_Water_Projection_mpi-esm1-2-hr_historical_0b7cdf4760.2023-11-22-131438', ...
              'GLOBE_Surface_Water_Projection_mpi-esm1-2-hr_ssp126_0b7cdf4760.2024-02-17-131701',     ...
              'GLOBE_Surface_Water_Projection_mpi-esm1-2-hr_ssp585_0b7cdf4760.2024-02-17-132324',     ...
              'GLOBE_Surface_Water_Projection_mri-esm2-0_historical_0b7cdf4760.2023-11-22-134259',    ...
              'GLOBE_Surface_Water_Projection_mri-esm2-0_ssp126_0b7cdf4760.2024-02-17-132852',        ...
              'GLOBE_Surface_Water_Projection_mri-esm2-0_ssp585_0b7cdf4760.2024-02-17-133437',        ...
              'GLOBE_Surface_Water_Projection_ukesm1-0-ll_historical_0b7cdf4760.2023-11-22-141413',   ...
              'GLOBE_Surface_Water_Projection_ukesm1-0-ll_ssp126_0b7cdf4760.2024-02-19-080920',       ...
              'GLOBE_Surface_Water_Projection_ukesm1-0-ll_ssp585_0b7cdf4760.2024-02-19-081526'};

for icase = 1 : 3%length(case_names)
    strs     = strsplit(case_names{icase},'_');
    model    = strs{5};
    scenario = strs{6};
    
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
    
    elmfiles = dir(['../outputs/' case_names{icase} '/run/*.elm.h0.*nc']);
    mosfiles = dir(['../outputs/' case_names{icase} '/run/*.mosart.h0.*nc']);

    elmfiles = elmfiles((yr_read-yr_start)*12+1:end);
    mosfiles = mosfiles((yr_read-yr_start)*12+1:end);
    assert(length(elmfiles) == length(mosfiles));
    assert(length(elmfiles) == (yr_end - yr_read + 1) * 12);
        
    if ~exist(['../data/projection_' model '_' scenario '.mat'],'file')
        
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
        
        save(['../data/projection_' model '_' scenario '.mat'],'fh2osfc','flooded','zwt','tws','rain','tsa','fsat','perch');
        
    elseif ~exist(['../data/projection_' model '_' scenario '_ice.mat'],'file')
    
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
        
        save(['../data/projection_' model '_' scenario '_ice.mat'],'qsnow','frost','fsno','snow','soilice','surfice');
        
    elseif ~exist(['../data/projection_' model '_' scenario '_flux.mat'],'file')
        
        for i = 1 : length(elmfiles)
            felm = fullfile(elmfiles(i).folder,elmfiles(i).name);
            fmos = fullfile(mosfiles(i).folder,mosfiles(i).name);
            
            disp([model ', ' scenario ': ' num2str(i) '/' num2str(length(elmfiles))]);
            
            if i == 1
                latent      = ncread(felm,'EFLX_LH_TOT');
            else
                latent(:,i) = ncread(felm,'EFLX_LH_TOT');

            end
        end
        
        save(['../data/projection_' model '_' scenario '_flux.mat'],'latent');
    
    end
    
    if ~exist(['../data/projection_' model '_' scenario '_energy.mat'],'file')
        
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
        
        save(['../data/projection_' model '_' scenario '_energy.mat'],'fsds','fsa','fsr','fire');
    
    end
    
end