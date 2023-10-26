clear;close all;clc;

% This script is to show the relationship between Surface water dynamics
% and Precipitation (Pr) and Temperature (Ta).

[e3sm_input, exportfig] = SetupEnvironment();

home = getenv('HOME');

days_of_month = [31; 28; 31; 30; 31; 30; 31; 31; 30; 31; 30; 31];

if strcmp(home,'/qfs/people/xudo627')
    read_gswp3 = 1;
    gswp3_dir = '/compyfs/inputdata/atm/datm7/atm_forcing.datm7.GSWP3.0.5d.v1.c170516/';
else
    read_gswp3 = 0;
end

if read_gswp3
    yr1 = 1993;
    yr2 = 2014;
    
    for yr = yr1 : yr2
        
        pr_mo = NaN(720,360,12);
        ta_mo = NaN(720,360,12);
        for mo = 1 : 12
            
            disp(['yr: ' num2str(yr) ', mo: ' num2str(mo)]);

            if mo < 10
            pr_fn = [gswp3_dir 'Precip3Hrly/clmforc.GSWP3.c2011.0.5x0.5.Prec.' num2str(yr) '-0' num2str(mo) '.nc'];
            ta_fn = [gswp3_dir 'TPHWL3Hrly/clmforc.GSWP3.c2011.0.5x0.5.TPQWL.' num2str(yr) '-0' num2str(mo) '.nc'];
            else
            pr_fn = [gswp3_dir 'Precip3Hrly/clmforc.GSWP3.c2011.0.5x0.5.Prec.' num2str(yr) '-'  num2str(mo) '.nc'];
            ta_fn = [gswp3_dir 'TPHWL3Hrly/clmforc.GSWP3.c2011.0.5x0.5.TPQWL.' num2str(yr) '-'  num2str(mo) '.nc'];
            end
            
            pr_mo(:,:,mo) = nanmean(ncread(pr_fn,'PRECmms')).*86400.*days_of_month(mo); % [mm/month]
            ta_mo(:,:,mo) = nanmean(ncread(ta_fn,'TBOT'),3); % [k]
        end

        pr_yr = nansum(pr,3); % [mm/yr]
        ta_mo = nanmean(ta_mo,3); % [k]

        out_fn = ['../data/GSWP3v1/GSWP3v1_' num2str(yr) '_Pr_Ta.mat'];
        save(out_fn,'pr_yr','ta_mo');

    end

end
