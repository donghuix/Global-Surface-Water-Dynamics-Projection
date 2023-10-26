clear;close all;clc;

[e3sm_input, exportfig] = SetupEnvironment();
domain_file = [e3sm_input 'share/domains/domain.lnd.r05_oEC60to30v3.190418.nc'];
re = 6.37122e6;% Earth radius

cmap = getPanoply_cMap('EO_aura_omi_formal');
cmap = create_nonlinear_cmap(cmap,0,0.2,0.2,3);

load('/Users/xudo627/Library/CloudStorage/OneDrive-PNNL/DATA/GIEMS/GIMES_half.mat');
load('../data/swf_cal.mat','swf_mon_cal');
xc   = ncread('../data/domain_lnd_GLOBE_1d.nc','xc');
yc   = ncread('../data/domain_lnd_GLOBE_1d.nc','yc');
xv   = ncread('../data/domain_lnd_GLOBE_1d.nc','xv');
yv   = ncread('../data/domain_lnd_GLOBE_1d.nc','yv');
area = ncread('../data/domain_lnd_GLOBE_1d.nc','area');

load('index_lnd.mat');
giems05 = NaN(length(xc),180);

for i = 1 : 180
    tmp = fliplr(GIEMS_half(:,:,i));
    giems05(:,i) = tmp(index_lnd);
end
load('LargeLakes.mat');
lakein = [];
for i = 1 : 20
    tmp = inpoly2([xc yc],[LargeLakes(i).X' LargeLakes(i).Y']);
    tmp = find(tmp == 1);
    lakein  = [lakein; tmp];
end

giems05(lakein,:)     = NaN;
swf_mon_cal(lakein,:) = NaN;
swf_mon_cal = swf_mon_cal(:,1:180);

giems05(isnan(swf_mon_cal)) = NaN;
swf_mon_cal(isnan(giems05)) = NaN;

continent = struct([]);
continent_code = {'af',    'ar',    'as',  'au',        'eu',    'gr',       'na',           'sa',           'si'     };
continent_name = {'Africa','Arctic','Asia','Austrialia','Europe','Greenland','North America','South America','Siberia'};

figure(1)
tmp = nansum(giems05.*area.*(re^2)./1e6,1);
subplot(3,3,1);
plot(zscore(tmp),'k-','LineWidth',2); hold on; grid on;
tmp = nanmean(reshape(tmp(:),[12,15]),1);
figure(2); set(gcf,'Position',[10 10 1200 600]);
subplot(3,3,1);
plot(1993:2007,zscore(tmp),'k-','LineWidth',2); hold on; grid on;

figure(1);
subplot(3,3,1);
tmp = nansum(swf_mon_cal.*area.*(re^2)./1e6,1);
plot(zscore(tmp),'r--','LineWidth',2);
tmp = nanmean(reshape(tmp(:),[12,15]),1);
add_title(gca,'Global');

figure(2);
subplot(3,3,1);
plot(1993:2007,zscore(tmp),'r--','LineWidth',2);
add_title(gca,'Global');
xlim([1993 2007]);

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
   tmp = nansum(giems05(continent(i).index,:).*area(continent(i).index).*(re^2)./1e6,1);
   subplot(3,3,k);
   plot(zscore(tmp),'k-','LineWidth',2); hold on; grid on;
   tmp = nanmean(reshape(tmp(:),[12,15]),1);
   figure(2);
   subplot(3,3,k);
   plot(1993:2007,zscore(tmp),'k-','LineWidth',2); hold on; grid on;
   
   figure(1);
   subplot(3,3,k);
   tmp = nansum(swf_mon_cal(continent(i).index,:).*area(continent(i).index).*(re^2)./1e6,1);
   plot(zscore(tmp),'r--','LineWidth',2);
   add_title(gca,continent_name{i});
   tmp = nanmean(reshape(tmp(:),[12,15]),1);
   figure(2);
   subplot(3,3,k);
   plot(1993:2007,zscore(tmp),'r--','LineWidth',2);
   add_title(gca,continent_name{i});
   xlim([1993 2007]);

   k = k + 1;

%    figure(2+i); set(gcf,'Position',[10 10 1200 500]);
%    subplot(1,2,1);
%    patch(xv(:,continent(i).index),yv(:,continent(i).index), ...
%          nanmean(giems05(continent(i).index,:),2),'LineStyle','none'); 
%    clim([0 0.2]); colormap(gca,cmap); hold on;
%    for j = 1 : length(S)
%        plot(S(j).X,S(j).Y,'k-','LineWidth',1); hold on;
%    end
% 
%    subplot(1,2,2);
%    patch(xv(:,continent(i).index),yv(:,continent(i).index), ...
%          nanmean(swf_mon_cal(continent(i).index,:),2),'LineStyle','none'); 
%    clim([0 0.2]); colormap(gca,cmap); hold on;
%    for j = 1 : length(S)
%        plot(S(j).X,S(j).Y,'k-','LineWidth',1); hold on;
%    end
%    R2(i) = estimate_evaluation_metric(nanmean(giems05(continent(i).index,:),2), ...
%                                       nanmean(swf_mon_cal(continent(i).index,:),2));
end

load coastlines.mat;
figure; set(gcf,'Position',[10 10 1000 1200]);
axs(1) = subplot(2,1,1);
patch(xv,yv,nanmean(giems05,2),'LineStyle','none'); hold on;
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
add_title(axs(1),'(a) GIEMS',20,'out');
add_title(axs(2),'(b) Simulation',20,'out');

load('colorblind_colormap.mat');
colorblind(5,:) = [];


% load('coastlines.mat');
% figure; set(gcf,'Position',[10 10 1200 500])
% subplot(1,2,1);
% patch(xv,yv,nanmean(giems05,2) - nanmin(giems05,[],2),'LineStyle','none'); clim([0 0.2]); colormap(gca,cmap); hold on;
% plot(coastlon,coastlat,'k-','LineWidth',1);
% ylim([-60 80]);
% 
% subplot(1,2,2);
% patch(xv,yv,nanmean(swf_mon_cal,2),'LineStyle','none'); clim([0 0.2]); colormap(gca,cmap); hold on;
% plot(coastlon,coastlat,'k-','LineWidth',1);
% ylim([-60 80]);
% 
% R2 = estimate_evaluation_metric(nanmean(giems05,2),nanmean(swf_mon_cal,2));

% giems025 = zeros(1440,720,180);
% data = load('/Users/xudo627/Library/CloudStorage/OneDrive-PNNL/DATA/GIEMS/wetland_global_extent_1993_2007_Papa_etal_2010_Prigent_etal_2012.dat');
% lat = data(:,2);
% lon = data(:,3);
% lon(lon > 180) = lon(lon > 180) - 360;
% giems = data(:,4:end);
% giems(giems == -99) = NaN;
% clear data;
% 
% lon025 = -180 + 0.25/2 : 0.25 : 180 - 0.25/2;
% lat025 = 90 - 0.25/2 : -0.25 : -90 + 0.25/2;
% [lon025, lat025] = meshgrid(lon025,lat025);
% lon025 = lon025'; lat025 = lat025';
% 
% for i = 1 : 180
%     tmp = data(:,i+3);
%     tmp(tmp == -99) = NaN; 
% end
% tmp = data(:,10); tmp(tmp == -99) = NaN; tmp(tmp < 0) = NaN; tmp = tmp ./ 773;
% 
% giems = griddata(lon,lat,tmp,lon025,lat025);
% 
% figure;
% [cb,~,~] = plot_globalspatial(lon025,lat025,giems,1,1);

