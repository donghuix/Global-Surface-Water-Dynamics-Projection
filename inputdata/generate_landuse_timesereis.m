clear;close all;clc;

addpath('/global/homes/d/donghui/Setup-E3SM-Mac/matlab-scripts-to-process-inputs');

dir = ['/global/cfs/cdirs/e3sm/inputdata/lnd/clm2/surfdata_map/'];
fnames = {'landuse.timeseries_0.5x0.5_hist_simyr1850-2015_c191004.nc',       ...
          'landuse.timeseries_0.5x0.5_SSP1_RCP26_simyr2015-2100_c220317.nc', ...
          'landuse.timeseries_0.5x0.5_SSP5_RCP85_simyr2015-2100_c220318.nc'};

out_netcdf_dir = '.';

lon_region = ncread('domain_lnd_GLOBE_1d.nc','xc');
lat_region = ncread('domain_lnd_GLOBE_1d.nc','yc');

%for i = 1 : length(fnames)
%    fname_in = [dir fnames{i}];
%    
%    if i == 1
%        usrdat_name = '0.5x0.5_HIST_simyr1850-2015_GLOBAL_1d';
%    elseif i == 2 
%        usrdat_name = '0.5x0.5_SSP1_RCP26_simyr2015-2100_GLOBAL_1d';
%    elseif i == 3
%        usrdat_name = '0.5x0.5_SSP5_RCP85_simyr2015-2100_GLOBAL_1d';
%    end
%    
%    fname_out = CreateELM_LandUse( lon_region, lat_region, ones(length(lon_region),1),   ...
%                                   fname_in, out_netcdf_dir, usrdat_name );
%end

fin  = 'landuse.timeseries_0.5x0.5_HIST_simyr1850-2015_GLOBAL_1d_c240229.nc';
fout = 'landuse.timeseries_0.5x0.5_HIST_simyr2015-2100_GLOBAL_1d_c240209.nc';
info = ncinfo(fin);
vars = {info.Variables.Name};

for i = 1 : length(vars)
    if strcmp(vars{i},'YEAR') || strcmp(vars{i},'time') || strcmp(vars{i},'input_pftdata_filename')
    elseif strcmp(vars{i},'PCT_NAT_PFT')
        varhist = ncread(fin,vars{i});
        varhist = repmat(nanmean(varhist(:,:,(1971 - 1850 + 1 : 2000 - 1850 + 1)),3),1,1,86);
        ncwrite(fout,vars{i},varhist);
        
    elseif strcmp(vars{i},'HARVEST_VH1') || strcmp(vars{i},'HARVEST_VH2') || strcmp(vars{i},'HARVEST_SH1') || ...
           strcmp(vars{i},'HARVEST_SH2') || strcmp(vars{i},'HARVEST_SH3') || strcmp(vars{i},'GRAZING')
        
        varhist = ncread(fin,vars{i});
        varhist = repmat(nanmean(varhist(:,(1971 - 1850 + 1 : 2000 - 1850 + 1)),2),1,86);
        ncwrite(fout,vars{i},varhist);
    else
        varhist = ncread(fin,vars{i});
        ncwrite(fout,vars{i},varhist);
    end
    
end

% 346.6991
fco2 = '/global/cfs/projectdirs/m3780/donghui/inputdata/CO2/fco2_datm_globalSSP5-8.5__simyr_2014-2501_CMIP6_c190506.nc';
fout = '/global/cfs/projectdirs/m3780/donghui/Global-Surface-Water-Dynamics-Projection/inputdata/fco2_datm_global_Hist_CO2__simyr_2014-2501_CMIP6_c190506.nc';

copyfile(fco2,fout);

CO2 = ncread(fout,'CO2');
CO2 = ones(size(CO2)).*364.6991;

ncwrite(fout,'CO2',CO2);
