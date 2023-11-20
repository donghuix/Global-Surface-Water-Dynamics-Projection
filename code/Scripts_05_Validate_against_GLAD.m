clear;close all;clc;

[e3sm_input, exportfig] = SetupEnvironment();

month_labels = {'01_Jan','02_Feb','03_Mar','04_Apr','05_May','06_Jun', ...
                '07_Jul','08_Aug','09_Sep','10_Oct','11_Nov','12_Dec'};
domain_file = [e3sm_input 'share/domains/domain.lnd.r05_oEC60to30v3.190418.nc'];
frac = ncread(domain_file,'frac'); frac = flipud(frac');
cmap = getPanoply_cMap('EO_aura_omi_formal');
cmap = create_nonlinear_cmap(cmap,0,0.2,0.2,3);
load('coastlines.mat');

% Read GLAD annual surface water dynamics at 0.5x0.5
load('../data/GLAD05/GLAD05_permenant_water.mat');
load('index_lnd.mat');
xc = ncread('../data/domain_lnd_GLOBE_1d.nc','xc');
yc = ncread('../data/domain_lnd_GLOBE_1d.nc','yc');
xv = ncread('../data/domain_lnd_GLOBE_1d.nc','xv');
yv = ncread('../data/domain_lnd_GLOBE_1d.nc','yv');
area = ncread('../data/domain_lnd_GLOBE_1d.nc','area');
re = 6.37122e6;% Earth radius
area = area.*(re^2) ./ 1e6; % square km

pw(frac < 1) = NaN;
gladseason = NaN(length(index_lnd),12);
gladannual = NaN(length(index_lnd),14);

for i = 2001 : 2014
    load(['../data/GLAD05/GLAD05_' num2str(i) '.mat']);
    disp(['Year-' num2str(i)]);

    GLAD05(frac < 1) = NaN;
    GLAD05 = GLAD05 - pw;
    GLAD05(GLAD05 <= 0) = 0;
    GLAD05 = fliplr(GLAD05');

    gladannual(:,i-2000)  = GLAD05(index_lnd);
end

for i = 1 : 12
    load(['../data/GLAD05/GLAD05_' month_labels{i} '.mat']);

    GLAD05(frac < 1) = NaN;
    GLAD05 = GLAD05 - pw;
    GLAD05(GLAD05 <= 0) = 0;
    GLAD05 = fliplr(GLAD05');

    gladseason(:,i)   = GLAD05(index_lnd);
end

pr1d = NaN(length(index_lnd),2014-2001+1);
ta1d = NaN(length(index_lnd),2014-2001+1);
for yr = 2001 : 2014
    load(['../data/GSWP3v1/GSWP3v1_' num2str(yr) '_Pr_Ta.mat']);
    tmp = pr_yr;
    pr_yr(1:360,:) = tmp(361:end,:);
    pr_yr(361:end,:) = tmp(1:360,:);
    pr1d(:,yr-2000) = pr_yr(index_lnd);

    tmp = ta_mo;
    ta_mo(1:360,:) = tmp(361:end,:);
    ta_mo(361:end,:) = tmp(1:360,:);
    ta1d(:,yr-2000) = ta_mo(index_lnd);
end

load('../data/swf_cal.mat','swf_yr_cal','swf_mon_cal');

load('LargeLakes.mat');
lakein = [];
for i = 1 : 20
    tmp = inpoly2([xc yc],[LargeLakes(i).X' LargeLakes(i).Y']);
    tmp = find(tmp == 1);
    lakein  = [lakein; tmp];
end

swf_yr_cal(lakein,:) = NaN; swf_yr_cal = swf_yr_cal(:,(2001 - 1993 + 1) : end);
gladannual(lakein,:) = NaN;
swf_yr_cal(nanmean(swf_mon_cal,2) >= 0.5,:) = NaN;
pr1d(lakein,:) = NaN; ta1d(lakein,:) = NaN;
pr1d(nanmean(swf_mon_cal,2) >= 0.5,:) = NaN;
ta1d(nanmean(swf_mon_cal,2) >= 0.5,:) = NaN;

swf_yr_cal(isnan(gladannual)) = NaN;
gladannual(isnan(swf_yr_cal)) = NaN;

swf_mon_cal(lakein,:) = NaN;
swf_mon_cal(nanmean(swf_mon_cal,2) >= 0.5,:) = NaN;

swf_mon_cal(isnan(gladseason)) = NaN;


swf_yr_cal  = swf_yr_cal  .* area;
swf_mon_cal = swf_mon_cal .* area;
pr1d        = pr1d .* area;
ta1d        = ta1d .* area;
gladannual  = gladannual./100 .* area;
gladseason  = gladseason./100 .* area;

continent = struct([]);
continent_code = {'af',    'ar',    'as',  'au',        'eu',    'gr',       'na',           'sa',           'si'     };
continent_name = {'Africa','Arctic','Asia','Austrialia','Europe','Greenland','North America','South America','Siberia'};
figure(1); set(gcf,'Position',[10 10 1200 600]);
subplot(3,3,1);
plot(2001:2014,zscore(nansum(gladannual,1)),'k-','LineWidth',2); hold on; grid on;
plot(2001:2014,zscore(nansum(swf_yr_cal,1)),'k--','LineWidth',2);
plot(2001:2014,zscore(nansum(pr1d,1)),'b--','LineWidth',2);
plot(2001:2014,zscore(nansum(ta1d,1)),'r--','LineWidth',2);
add_title(gca,'Global');
xlim([2001 2014]);

y = zscore(nansum(gladannual,1));
x1= zscore(nansum(swf_yr_cal,1));
x2= zscore(nansum(pr1d,1));
x3= zscore(nansum(ta1d,1));

y = detrend(y);
x1= detrend(x1);
x2= detrend(x2);
x3= detrend(x3);

R = corrcoef(y,x1);
rho(1,1) = R(1,2);
R = corrcoef(x1,x2);
rho(1,2) = R(1,2);
R = corrcoef(x1,x3);
rho(1,3) = R(1,2);
R = corrcoef(x2,x3);
rho(1,4) = R(1,2);


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
   plot(2001:2014,zscore(nansum(gladannual(continent(i).index,:),1)),'k-','LineWidth',2); hold on; grid on;
   plot(2001:2014,zscore(nansum(swf_yr_cal(continent(i).index,:),1)),'k--','LineWidth',2);
   plot(2001:2014,zscore(nansum(pr1d(continent(i).index,:),1)),'b--','LineWidth',2);
   plot(2001:2014,zscore(nansum(ta1d(continent(i).index,:),1)),'r--','LineWidth',2);

   y = zscore(nansum(gladannual(continent(i).index,:),1));
   x1= zscore(nansum(swf_yr_cal(continent(i).index,:),1));
   x2= zscore(nansum(pr1d(continent(i).index,:),1));
   x3= zscore(nansum(ta1d(continent(i).index,:),1));

   y = detrend(y);
   x1= detrend(x1);
   x2= detrend(x2);
   x3= detrend(x3);
   
   R = corrcoef(y,x1);
   rho(k,1) = R(1,2);
   R = corrcoef(x1,x2);
   rho(k,2) = R(1,2);
   R = corrcoef(x1,x3);
   rho(k,3) = R(1,2);
   R = corrcoef(x2,x3);
   rho(k,4) = R(1,2);

   add_title(gca,continent_name{i});
   xlim([2001 2014]);

   figure(2);
   subplot(3,3,i);
   plot(nansum(gladseason(continent(i).index,:),1), 'k-','LineWidth',2); hold on; grid on;
   plot(nansum(swf_mon_cal(continent(i).index,:),1),'r--','LineWidth',2);
   add_title(gca,continent_name{i});
   k = k + 1;
end

% i = 2;
% annual = gladannual(continent(i).index,1);
% annual2 = nanmean(swf_yr_cal(continent(i).index,:),2);
% month  = nanmean(gladseason(continent(i).index,:),2);
% month2  = nanmean(swf_mon_cal(continent(i).index,:),2);
% 
% tmpxv  = xv(:,continent(i).index);
% tmpyv  = yv(:,continent(i).index);
% 
% figure;
% patch(tmpxv,tmpyv,annual,'LineStyle','none'); colorbar;
% figure;
% patch(tmpxv,tmpyv,month,'LineStyle','none');colorbar;
% 
% figure;
% patch(tmpxv,tmpyv,annual2,'LineStyle','none'); colorbar; 
% figure;
% patch(tmpxv,tmpyv,month2,'LineStyle','none');colorbar;
% 
% 
% figure;
% plot(nansum(gladannual(continent(i).index,:),1),'k-','LineWidth',2); hold on; grid on;
% plot(nansum(swf_yr_cal(continent(i).index,:),1),'r--','LineWidth',2);
% 
% figure;
% plot(nansum(gladseason(continent(i).index,:),1), 'k-','LineWidth',2); hold on; grid on;
% plot(nansum(swf_mon_cal(continent(i).index,:),1),'r--','LineWidth',2);