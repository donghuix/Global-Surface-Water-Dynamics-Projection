clear;close all;clc;

addpath('/global/homes/d/donghui/Setup-E3SM-Mac/matlab-scripts-to-process-inputs');

dir = ['/global/cfs/cdirs/e3sm/inputdata/lnd/clm2/surfdata_map/'];
fnames = {'landuse.timeseries_0.5x0.5_SSP1_RCP26_simyr2015-2100_c220317.nc', ...
          'landuse.timeseries_0.5x0.5_SSP5_RCP85_simyr2015-2100_c220318.nc'};

out_netcdf_dir = '.';

lon_region = ncread('domain_lnd_GLOBE_1d.nc','xc');
lat_region = ncread('domain_lnd_GLOBE_1d.nc','yc');

for i = 1 : length(fnames)
    fname_in = [dir fnames{i}];
    
    if i == 1 
        usrdat_name = '0.5x0.5_SSP1_RCP26_simyr2015-2100_GLOBAL_1d';
    elseif i == 2
        usrdat_name = '0.5x0.5_SSP5_RCP85_simyr2015-2100_GLOBAL_1d';
    end
    
    fname_out = CreateELM_LandUse( lon_region, lat_region, ones(length(lon_region),1),   ...
                                   fname_in, out_netcdf_dir, usrdat_name );
end

