clear;close all;clc;

[e3sm_input, exportfig] = SetupEnvironment();

domain_file = [e3sm_input 'share/domains/domain.lnd.r05_oEC60to30v3.190418.nc'];
frac = ncread(domain_file,'frac'); %frac = flipud(frac');

cmap = getPanoply_cMap('EO_aura_omi_formal');
cmap = create_nonlinear_cmap(cmap,0,0.2,0.2,3);

load('index_lnd.mat');
xc = ncread('../data/domain_lnd_GLOBE_1d.nc','xc');
yc = ncread('../data/domain_lnd_GLOBE_1d.nc','yc');
xv = ncread('../data/domain_lnd_GLOBE_1d.nc','xv');
yv = ncread('../data/domain_lnd_GLOBE_1d.nc','yv');
area = ncread('../data/domain_lnd_GLOBE_1d.nc','area');
re = 6.37122e6;% Earth radius
area = area.*(re^2) ./ 1e6; % square km

lon  = ncread('../data/WAD2M_wetlands_2000-2020_05deg_Ver2.0.nc','lon');
lat  = ncread('../data/WAD2M_wetlands_2000-2020_05deg_Ver2.0.nc','lat');
Fw05 = ncread('../data/WAD2M_wetlands_2000-2020_05deg_Ver2.0.nc','Fw'); 
Fw05 = Fw05(:,:,(2001-2000)*12+1 : (2014-2000+1)*12);

figure;
imagesc(nanmean(Fw05,3));
figure;
imagesc(frac);
nt = size(Fw05,3);

Fw1d = NaN(length(index_lnd),nt);
for i = 1 : nt
    tmp = Fw05(:,:,i);
    Fw1d(:,i) = tmp(index_lnd);
end

load('../data/swf_cal.mat','swf_yr_cal','swf_mon_cal');

load('LargeLakes.mat');
lakein = [];
for i = 1 : 20
    tmp = inpoly2([xc yc],[LargeLakes(i).X' LargeLakes(i).Y']);
    tmp = find(tmp == 1);
    lakein  = [lakein; tmp];
end

Fw1d(lakein,:)     = NaN;
swf_mon_cal(lakein,:) = NaN;
swf_mon_cal = swf_mon_cal(:,(2001-1993)*12+1:(2014-1993+1)*12);

continent = struct([]);
continent_code = {'af',    'ar',    'as',  'au',        'eu',    'gr',       'na',           'sa',           'si'     };
continent_name = {'Africa','Arctic','Asia','Austrialia','Europe','Greenland','North America','South America','Siberia'};

xt = 2001 : 2014;
figure(1)
tmp = nansum(Fw1d.*area.*(re^2)./1e6,1);
subplot(3,3,1);
plot(zscore(tmp),'k-','LineWidth',2); hold on; grid on;
tmp = nanmean(reshape(tmp(:),[12,14]),1);
figure(2); set(gcf,'Position',[10 10 1200 600]);
subplot(3,3,1);
plot(xt,zscore(tmp),'k-','LineWidth',2); hold on; grid on;

figure(1);
subplot(3,3,1);
tmp = nansum(swf_mon_cal.*area.*(re^2)./1e6,1);
plot(zscore(tmp),'r--','LineWidth',2);
tmp = nanmean(reshape(tmp(:),[12,14]),1);
add_title(gca,'Global');

figure(2);
subplot(3,3,1);
plot(xt,zscore(tmp),'r--','LineWidth',2);
add_title(gca,'Global');
xlim([xt(1) xt(end)]);

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
      continent(i).index = [continent(i).index; tmp];
   end

   figure(1)
   tmp = nansum(Fw1d(continent(i).index,:).*area(continent(i).index).*(re^2)./1e6,1);
   subplot(3,3,k);
   plot(zscore(tmp),'k-','LineWidth',2); hold on; grid on;
   tmp = nanmean(reshape(tmp(:),[12,14]),1);
   figure(2);
   subplot(3,3,k);
   plot(xt,zscore(tmp),'k-','LineWidth',2); hold on; grid on;
   
   figure(1);
   subplot(3,3,k);
   tmp = nansum(swf_mon_cal(continent(i).index,:).*area(continent(i).index).*(re^2)./1e6,1);
   plot(zscore(tmp),'r--','LineWidth',2);
   add_title(gca,continent_name{i});
   tmp = nanmean(reshape(tmp(:),[12,14]),1);
   figure(2);
   subplot(3,3,k);
   plot(xt,zscore(tmp),'r--','LineWidth',2);
   add_title(gca,continent_name{i});
   xlim([xt(1) xt(end)]);

   k = k + 1;
end

load coastlines.mat;
figure; set(gcf,'Position',[10 10 1000 1200]);
axs(1) = subplot(2,1,1);
patch(xv,yv,nanmean(Fw1d,2),'LineStyle','none'); hold on;
clim([0 0.3]); colormap(gca,cmap); ylim([-60 80]); xlim([-180 180]);
plot(coastlon,coastlat,'k-','LineWidth',1.5);
axs(2) = subplot(2,1,2);
patch(xv,yv,nanmean(swf_mon_cal,2),'LineStyle','none'); hold on;
clim([0 0.3]); colormap(gca,cmap); ylim([-60 80]); xlim([-180 180]);
plot(coastlon,coastlat,'k-','LineWidth',1.5);
cb = colorbar('south');

axs(2).Position(2) = axs(2).Position(2) + 0.075;
axs(1).Position(1) = axs(1).Position(1) - 0.1;
axs(2).Position(1) = axs(2).Position(1) - 0.1;
axs(1).Position(3) = axs(1).Position(3) - 0.1;
axs(2).Position(3) = axs(2).Position(3) - 0.1;

cb.Position(1) = axs(2).Position(1);
cb.Position(2) = axs(2).Position(2)-0.05;
cb.Position(3) = axs(2).Position(3);
cb.Position(4) = 0.02;
cb.AxisLocation = 'out'; cb.FontSize = 16;
add_title(axs(1),'(a) WAD2M',20,'out');
add_title(axs(2),'(b) Simulation',20,'out');

load('colorblind_colormap.mat');
colorblind(5,:) = [];

[R2,RMSE,NSE,PBIAS,MSE,NSE1] = estimate_evaluation_metric(nanmean(Fw1d,2),nanmean(swf_mon_cal,2));

