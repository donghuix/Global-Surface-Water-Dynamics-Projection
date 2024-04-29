clear;close all;clc;

[e3sm_input, exportfig] = SetupEnvironment();

model = 'gfdl-esm4';
mid_start = 2041;
mid_end   = 2070;
end_start = 2071;
end_end   = 2100;
month_labels = {'01_Jan','02_Feb','03_Mar','04_Apr','05_May','06_Jun', ...
                '07_Jul','08_Aug','09_Sep','10_Oct','11_Nov','12_Dec'};

load('index_lnd.mat');
xv = ncread('../data/domain_lnd_GLOBE_1d.nc','xv');
yv = ncread('../data/domain_lnd_GLOBE_1d.nc','yv');
xc = ncread('../data/domain_lnd_GLOBE_1d.nc','xc');
yc = ncread('../data/domain_lnd_GLOBE_1d.nc','yc');
area = ncread('../data/domain_lnd_GLOBE_1d.nc','area');
domain_file = [e3sm_input 'share/domains/domain.lnd.r05_oEC60to30v3.190418.nc'];
frac = ncread(domain_file,'frac'); frac = flipud(frac');
load('../data/GLAD05/GLAD05_permenant_water.mat');
pw(frac < 1) = NaN;
gladseason = NaN(length(index_lnd),12);
gladannual = NaN(720,360,12);
for i = 1 : 12
    load(['../data/GLAD05/GLAD05_' month_labels{i} '.mat']);
    GLAD05(frac < 1) = NaN;
    GLAD05 = GLAD05 - pw;
    GLAD05(GLAD05 <= 0) = 0;
    GLAD05 = fliplr(GLAD05');

    gladseason(:,i)   = GLAD05(index_lnd);
    gladannual(:,:,i) = GLAD05;
    patch(xv,yv,GLAD05(index_lnd)./100,'LineStyle','none'); clim([0 0.2]); hold on; grid on;
    ylim([-60 80]); set(gca,'XTick',[],'YTick',[]);
end
gladannual = nanmean(gladannual,3);

re = 6.37122e6;             % Earth radius
area = area.*(re^2) ./ 1e6; % square km
[lakein,lakein2d] = getLakeIndex(e3sm_input);

hist  = load(['../projection/projection_BC_' model '_historical.mat'],'flooded','fh2osfc','rain');
futBD = load(['../projection/projection_BC_' model '_ssp585.mat'],'flooded','fh2osfc','rain');
hisBD = load(['../projection/projection_histBD_' model '_ssp585.mat'],'flooded','fh2osfc');
hisLU = load(['../projection/projection_histLU_' model '_ssp585.mat'],'flooded','fh2osfc');
hisCO = load(['../projection/projection_histCO2_' model '_ssp585.mat'],'flooded','fh2osfc');

hist.swf      = hist.flooded(:,1:360) + hist.fh2osfc(:,1:360);
hist.rain     = hist.rain(:,1:360);
futBD.end.rain= futBD.rain(:,(end_start-2015)*12+1:(end_end-2015+1)*12);
futBD.end.swf = futBD.flooded(:,(end_start-2015)*12+1:(end_end-2015+1)*12) + ...
                futBD.fh2osfc(:,(end_start-2015)*12+1:(end_end-2015+1)*12);
hisBD.end.swf = hisBD.flooded(:,(end_start-2015)*12+1:(end_end-2015+1)*12) + ...
                hisBD.fh2osfc(:,(end_start-2015)*12+1:(end_end-2015+1)*12);
hisLU.end.swf = hisLU.flooded(:,(end_start-2015)*12+1:(end_end-2015+1)*12) + ...
                hisLU.fh2osfc(:,(end_start-2015)*12+1:(end_end-2015+1)*12);
hisCO.end.swf = hisCO.flooded(:,(end_start-2015)*12+1:(end_end-2015+1)*12) + ...
                hisCO.fh2osfc(:,(end_start-2015)*12+1:(end_end-2015+1)*12);

hist.swf(lakein,:)      = NaN; hist.rain(lakein,:)  = NaN;
futBD.end.swf(lakein,:) = NaN; futBD.rain(lakein,:) = NaN;
hisBD.end.swf(lakein,:) = NaN;
hisLU.end.swf(lakein,:) = NaN;
hisCO.end.swf(lakein,:) = NaN;

cmap = getPanoply_cMap('UKM_hadcrut_10');
cmap(7,3) = 1;
cmap(12,:)= [];
cmap(6,:)= [];
cmap = flipud(cmap);

N = 2;
cmap_new = NaN(15*2*N+1,3);
cmap_new(1:5*N,:)  = repmat(cmap(1,:),5*N,1);
cmap_new(5*N+1:9*N,:) = repmat(cmap(2,:),4*N,1);
cmap_new(9*N+1:12*N,:) = repmat(cmap(3,:),3*N,1);
cmap_new(12*N+1:14*N,:) = repmat(cmap(4,:),2*N,1);
cmap_new(14*N+1:15*N,:) = repmat(cmap(5,:),N,1);
cmap_new(15*N+1,:)    = repmat(cmap(6,:),1,1);
cmap_new(15*N+2:16*N+1,:) = repmat(cmap(7,:),N,1);
cmap_new(16*N+2:18*N+1,:) = repmat(cmap(8,:),2*N,1);
cmap_new(18*N+2:21*N+1,:) = repmat(cmap(9,:),3*N,1);
cmap_new(21*N+2:25*N+1,:) = repmat(cmap(10,:),4*N,1);
cmap_new(25*N+2:30*N+1,:) = repmat(cmap(11,:),5*N,1);

load('coastlines.mat');

figure; set(gcf,'Position',[10 10 1200 600])
patch(xv,yv,nanmean(futBD.end.swf - hisCO.end.swf,2).*area,'LineStyle','none'); 
clim([-10 10]); colormap(cmap_new);hold on; grid on;
ylim([-60 80]); xlim([-180 180]);%set(gca,'XTick',[],'YTick',[]); 
plot(coastlon,coastlat,'k-','LineWidth',1);
ax1 = gca;
ax1.Position(1) = ax1.Position(1) - 0.1;
ax1.Position(3) = ax1.Position(3) - 0.05;
t = add_title(ax1,'Impacts of CO_{2} on Surface Water under SSP585 [km^2]',20,'out');
pos1 = ax1.Position;
cb = colorbar('south');
cb.Position(1) = pos1(1);
cb.Position(2) = pos1(2) - 0.1;
cb.Position(3) = pos1(3);
cb.FontSize = 15;

CO2d = NaN(720,360);
CO2d(index_lnd) = nanmean(futBD.end.swf - hisCO.end.swf,2).*area;
CO2d(isnan(gladannual)) = NaN;
ax2 = axes('Position',[pos1(1) + pos1(3) + 0.02 pos1(2) 0.15 ax1.Position(4)]);
plot(ax2,movmean(nansum(CO2d,1),5),-89.75:0.5:89.75,'b-', 'LineWidth',2);hold on; grid on;
ylim([-60 80]);


figure; set(gcf,'Position',[10 10 1200 600])
patch(xv,yv,nanmean(futBD.end.swf - hisLU.end.swf,2).*area,'LineStyle','none'); 
clim([-10 10]); colormap(cmap_new);hold on; grid on;
ylim([-60 80]); xlim([-180 180]);%set(gca,'XTick',[],'YTick',[]); 
plot(coastlon,coastlat,'k-','LineWidth',1);
ax1 = gca;
ax1.Position(1) = ax1.Position(1) - 0.1;
ax1.Position(3) = ax1.Position(3) - 0.05;
t = add_title(ax1,'Impacts of LUC on Surface Water under SSP585 [km^2]',20,'out');
pos1 = ax1.Position;
cb = colorbar('south');
cb.Position(1) = pos1(1);
cb.Position(2) = pos1(2) - 0.1;
cb.Position(3) = pos1(3);
cb.FontSize = 15;

LU2d = NaN(720,360);
LU2d(index_lnd) = nanmean(futBD.end.swf - hisLU.end.swf,2).*area;
LU2d(isnan(gladannual)) = NaN;
ax2 = axes('Position',[pos1(1) + pos1(3) + 0.02 pos1(2) 0.15 ax1.Position(4)]);
plot(ax2,movmean(nansum(LU2d,1),5),-89.75:0.5:89.75,'b-', 'LineWidth',2);hold on; grid on;
ylim([-60 80]);


figure; set(gcf,'Position',[10 10 1200 600])
patch(xv,yv,nanmean(futBD.end.swf - hisBD.end.swf,2).*area,'LineStyle','none'); 
clim([-10 10]); colormap(cmap_new);hold on; grid on;
ylim([-60 80]); xlim([-180 180]);%set(gca,'XTick',[],'YTick',[]); 
plot(coastlon,coastlat,'k-','LineWidth',1);
ax1 = gca;
ax1.Position(1) = ax1.Position(1) - 0.1;
ax1.Position(3) = ax1.Position(3) - 0.05;
t = add_title(ax1,'Impacts of BD on Surface Water under SSP585 [km^2]',20,'out');
pos1 = ax1.Position;
cb = colorbar('south');
cb.Position(1) = pos1(1);
cb.Position(2) = pos1(2) - 0.1;
cb.Position(3) = pos1(3);
cb.FontSize = 15;

BD2d = NaN(720,360);
BD2d(index_lnd) = nanmean(futBD.end.swf - hisBD.end.swf,2).*area;
BD2d(isnan(gladannual)) = NaN;
CC2d = NaN(720,360);
CC2d(index_lnd) = nanmean(futBD.end.swf - hist.swf,2).*area;
CC2d(isnan(gladannual)) = NaN;
ax2 = axes('Position',[pos1(1) + pos1(3) + 0.02 pos1(2) 0.15 ax1.Position(4)]);
plot(ax2,movmean(nansum(BD2d,1),5),-89.75:0.5:89.75,'b-', 'LineWidth',2);hold on; grid on;
ylim([-60 80]);

nansum(nanmean(futBD.end.swf - hist.swf,2).*area);
nansum(nanmean(futBD.end.swf - hisBD.end.swf,2).*area)

figure; set(gcf,'Position',[10 10 1200 600])
patch(xv,yv,nanmean(futBD.end.swf - hist.swf,2).*area,'LineStyle','none'); 
clim([-30.5 30.5]); colormap(cmap_new);hold on; grid on;
ylim([-60 80]); xlim([-180 180]);%set(gca,'XTick',[],'YTick',[]); 
plot(coastlon,coastlat,'k-','LineWidth',1);
ax1 = gca;
ax1.Position(1) = ax1.Position(1) - 0.1;
ax1.Position(3) = ax1.Position(3) - 0.05;
t = add_title(ax1,'Impacts of BD on Surface Water under SSP585',20,'out');
pos1 = ax1.Position;
ax2 = axes('Position',[pos1(1) + pos1(3) + 0.02 pos1(2) 0.15 ax1.Position(4)]);
plot(ax2,movmean(nansum(CC2d,1),5),-89.75:0.5:89.75,'b-', 'LineWidth',2);hold on; grid on;
plot(ax2,movmean(nansum(BD2d,1),5),-89.75:0.5:89.75,'r-', 'LineWidth',2);hold on; grid on;
plot(ax2,movmean(nansum(LU2d,1),5),-89.75:0.5:89.75,'g-', 'LineWidth',2);hold on; grid on;
plot(ax2,movmean(nansum(CO2d,1),5),-89.75:0.5:89.75,'k-', 'LineWidth',2);hold on; grid on;
legend('Total change','BD impacts','LU impacts','CO2 impacts','FontSize',15,'FontWeight','bold');
ylim([-60 80]);
pos1 = ax1.Position;
cb = colorbar(ax1,'south');
cb.Position(1) = pos1(1);
cb.Position(2) = pos1(2) - 0.1;
cb.Position(3) = pos1(3);
cb.FontSize = 15;
cb.Ticks = [-30.5 -20.5 -12.5 -6.5 -2.5 0 2.5 6.5 12.5 20.5 30.5];
cb.TickLabels = {'-30', '-20','-12','-6','-2','0','2','6','12','20','30'};
