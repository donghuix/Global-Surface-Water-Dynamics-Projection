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
            
            pr_mo(:,:,mo) = nanmean(ncread(pr_fn,'PRECTmms'),3).*86400.*days_of_month(mo); % [mm/month] h
            ta_mo(:,:,mo) = nanmean(ncread(ta_fn,'TBOT'),3); % [k]
        end

        pr_yr = nansum(pr_mo,3); % [mm/yr]
        ta_mo = nanmean(ta_mo,3); % [k]

        out_fn = ['../data/GSWP3v1/GSWP3v1_' num2str(yr) '_Pr_Ta.mat'];
        save(out_fn,'pr_yr','ta_mo');

    end

end

load('index_lnd.mat');
xc = ncread('../data/domain_lnd_GLOBE_1d.nc','xc');
yc = ncread('../data/domain_lnd_GLOBE_1d.nc','yc');
xv = ncread('../data/domain_lnd_GLOBE_1d.nc','xv');
yv = ncread('../data/domain_lnd_GLOBE_1d.nc','yv');
area = ncread('../data/domain_lnd_GLOBE_1d.nc','area');
re = 6.37122e6;% Earth radius
area = area.*(re^2) ./ 1e6; % square km

yr1 = 1993;
yr2 = 2014;

pr1d = NaN(length(index_lnd),yr2 - yr1 + 1);
ta1d = NaN(length(index_lnd),yr2 - yr1 + 1);
for yr = yr1 : yr2
    load(['../data/GSWP3v1/GSWP3v1_' num2str(yr) '_Pr_Ta.mat']);
    tmp = pr_yr;
    pr_yr(1:360,:) = tmp(361:end,:);
    pr_yr(361:end,:) = tmp(1:360,:);
    pr1d(:,yr-1993+1) = pr_yr(index_lnd);

    tmp = ta_mo;
    ta_mo(1:360,:) = tmp(361:end,:);
    ta_mo(361:end,:) = tmp(1:360,:);
    ta1d(:,yr-1993+1) = ta_mo(index_lnd);
end

load('../data/swf_cal.mat','swf_yr_cal','swf_mon_cal');
load('LargeLakes.mat');
lakein = [];
for i = 1 : 20
    tmp = inpoly2([xc yc],[LargeLakes(i).X' LargeLakes(i).Y']);
    tmp = find(tmp == 1);
    lakein  = [lakein; tmp];
end

swf_yr_cal(lakein,:) = NaN; 
swf_yr_cal(nanmean(swf_mon_cal,2) >= 0.5,:) = NaN;
pr1d(lakein,:) = NaN; ta1d(lakein,:) = NaN;
pr1d(nanmean(swf_mon_cal,2) >= 0.5,:) = NaN;
ta1d(nanmean(swf_mon_cal,2) >= 0.5,:) = NaN;

swf_yr_cal  = swf_yr_cal  .* area;
pr1d        = pr1d .* area;
ta1d        = ta1d .* area;

continent = struct([]);
continent_code = {'af',    'ar',    'as',  'au',        'eu',    'gr',       'na',           'sa',           'si'     };
continent_name = {'Africa','Arctic','Asia','Austrialia','Europe','Greenland','North America','South America','Siberia'};
figure(1); set(gcf,'Position',[10 10 1200 600]);
subplot(3,3,1);
plot(yr1 : yr2,zscore(nansum(swf_yr_cal,1)),'k--','LineWidth',2); hold on; grid on;
plot(yr1 : yr2,zscore(nansum(pr1d,1)),'b--','LineWidth',2);
plot(yr1 : yr2,zscore(nanmean(ta1d,1)),'r--','LineWidth',2);
add_title(gca,'Global');
xlim([yr1 yr2]);

x1= zscore(nansum(swf_yr_cal,1));
x2= zscore(nansum(pr1d,1));
x3= zscore(nanmean(ta1d,1));

% x1= detrend(x1);
% x2= detrend(x2);
% x3= detrend(x3);

R = corrcoef(x1,x2);
rho(1,1) = R(1,2);
R = corrcoef(x1,x3);
rho(1,2) = R(1,2);
R = corrcoef(x2,x3);
rho(1,3) = R(1,2);


k = 2;
for i = [1 2 3 4 5 7 8 9]
   code = continent_code{i};
   continent(i).code  = code;
   continent(i).name  = continent_name{i};
   continent(i).index = []; 
   continent(i).swf   = NaN(20,1);
   S = shaperead(['../data/HydroBASINS/hybas_' code '_lev01-06_v1c/hybas_' code '_lev01_v1c.shp']);
   for j = 1 : length(S)
       tmp = inpoly2([xc(:) yc(:)],[S(j).X' S(j).Y']);
       tmp = find(tmp == 1);
       %plot(S(j).X,S(j).Y,'-','LineWidth',2);
       continent(i).index = [continent(i).index; tmp];
   end
   
   figure(1);
   subplot(3,3,k);
   plot(yr1 : yr2,zscore(nansum(swf_yr_cal(continent(i).index,:),1)),'k--','LineWidth',2); hold on; grid on;
   plot(yr1 : yr2,zscore(nansum(pr1d(continent(i).index,:),1)),'b--','LineWidth',2);
   plot(yr1 : yr2,zscore(nanmean(ta1d(continent(i).index,:),1)),'r--','LineWidth',2);

   x1= zscore(nansum(swf_yr_cal(continent(i).index,:),1));
   x2= zscore(nansum(pr1d(continent(i).index,:),1));
   x3= zscore(nanmean(ta1d(continent(i).index,:),1));

%    x1= detrend(x1);
%    x2= detrend(x2);
%    x3= detrend(x3);
   
   R = corrcoef(x1,x2);
   rho(k,1) = R(1,2);
   R = corrcoef(x1,x3);
   rho(k,2) = R(1,2);
   R = corrcoef(x2,x3);
   rho(k,3) = R(1,2);

   add_title(gca,continent_name{i});
   xlim([yr1 yr2]);

   k = k + 1;
end

rho1d = NaN(length(xc),4);
xtrain = NaN(length(xc),22,2);
ytrain = NaN(length(xc),22);
for i = 1 : length(xc)
    
    disp(i);

    x1= zscore(swf_yr_cal(i,:));
    x2= zscore(pr1d(i,:));
    x3= zscore(ta1d(i,:));

    x1 = detrend(x1); x2 = detrend(x2); x3 = detrend(x3);

    R = corrcoef(x1,x2);
    rho1d(i,1) = R(1,2);
    R = corrcoef(x1,x3);
    rho1d(i,2) = R(1,2);
    R = corrcoef(x2,x3);
    rho1d(i,4) = R(1,2);

%     [b,bint,r,rint,stats] = regress(x1',[ones(size(x1')) x2' x3' x2'.*x3']);
%     rho1d(i,3) = sqrt(stats(1));
%     
    xtrain(i,:,1) = x2;
    xtrain(i,:,2) = x3;
    ytrain(i,:)   = x1;
end

save('train.mat','xtrain','ytrain');

figure; set(gcf,'Position',[10 10 800 1000])
axs(1) = subplot(3,1,1);
patch(xv,yv,rho1d(:,1),'LineStyle','none'); clim([-1.1 1.1]); colormap(blue2red(11));
axs(2) = subplot(3,1,2);
patch(xv,yv,rho1d(:,2),'LineStyle','none'); clim([-1.1 1.1]); colormap(blue2red(11));
axs(3) = subplot(3,1,3);
patch(xv,yv,rho1d(:,3),'LineStyle','none'); clim([-1.1 1.1]); colormap(blue2red(11));
cb = colorbar('east');

axs(1).Position(3) = axs(1).Position(3) - 0.1;
axs(2).Position(3) = axs(2).Position(3) - 0.1;
axs(3).Position(3) = axs(3).Position(3) - 0.1;
axs(1).Position(4) = axs(1).Position(4) + 0.05;
axs(2).Position(4) = axs(2).Position(4) + 0.05;
axs(3).Position(4) = axs(3).Position(4) + 0.05;
cb.Position(1) = axs(3).Position(1) + axs(3).Position(3) + 0.02;
cb.Position(2) = axs(3).Position(2);
cb.Position(3) = 0.02;
cb.Position(4) = axs(1).Position(2) + axs(1).Position(4) - cb.Position(2);
cb.AxisLocation = 'out'; cb.FontSize = 15;

rf = load('rf.dat');
rf(rf < -90) = NaN;
figure; set(gcf,'Position',[10 10 800 1000])
axs(1) = subplot(3,1,1);
patch(xv,yv,rf(:,1),'LineStyle','none'); clim([-1.1 1.1]); colormap(gca,blue2red(11));
axs(2) = subplot(3,1,2);
patch(xv,yv,rf(:,2),'LineStyle','none'); clim([0 1]); colormap(gca,jet);
axs(3) = subplot(3,1,3);
patch(xv,yv,rf(:,3),'LineStyle','none'); clim([0 1]); colormap(gca,jet);