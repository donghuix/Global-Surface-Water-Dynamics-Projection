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

re = 6.37122e6;% Earth radius
area = area.*(re^2) ./ 1e6; % square km
[lakein,lakein2d] = getLakeIndex(e3sm_input);

historical = load(['../projection/projection_cal12_' model '_historical.mat'],'flooded','fh2osfc','perch','zwt','tsa','rain','fsat');
ssp126     = load(['../projection/projection_cal12_' model '_ssp126.mat'],'flooded','fh2osfc','perch','zwt','tsa','rain','fsat');
ssp585     = load(['../projection/projection_cal12_' model '_ssp585.mat'],'flooded','fh2osfc','perch','zwt','tsa','rain','fsat');
flux.historical = load(['../projection/projection_cal12_' model '_historical_flux.mat'],'fgr','surftsoi','soiltemp');
flux.ssp126     = load(['../projection/projection_cal12_' model '_ssp126_flux.mat'],'fgr','surftsoi','soiltemp');
flux.ssp585     = load(['../projection/projection_cal12_' model '_ssp585_flux.mat'],'fgr','surftsoi','soiltemp');
ice.historical = load(['../projection/projection_cal12_' model '_historical_ice.mat'],'soilice','surfice','snow','fsno');
ice.ssp126     = load(['../projection/projection_cal12_' model '_ssp126_ice.mat'],'soilice','surfice','snow','fsno');
ice.ssp585     = load(['../projection/projection_cal12_' model '_ssp585_ice.mat'],'soilice','surfice','snow','fsno');

historical.flooded = historical.flooded(:,1:360);
historical.fh2osfc = historical.fh2osfc(:,1:360);
historical.swf     = historical.flooded + historical.fh2osfc;
historical.fsat    = historical.fsat(:,1:360);
historical.perch   = historical.perch(:,1:360);
historical.zwt     = historical.zwt(:,1:360);
historical.tsa     = historical.tsa(:,1:360)- 273.15;
historical.rain    = historical.rain(:,1:360);
flux.historical.fgr    = flux.historical.fgr(:,1:360);
flux.historical.surftsoi = flux.historical.surftsoi(:,1:360);
flux.historical.soiltemp = flux.historical.soiltemp(:,1:360);
ice.historical.soilice = ice.historical.soilice(:,1:360);
ice.historical.surfice = ice.historical.surfice(:,1:360);
ice.historical.snow    = ice.historical.snow(:,1:360);
ice.historical.fsno    = ice.historical.fsno(:,1:360);


ssp126.mid.flooded = ssp126.flooded(:,(mid_start-2015)*12+1:(mid_end-2015+1)*12);
ssp126.mid.fh2osfc = ssp126.fh2osfc(:,(mid_start-2015)*12+1:(mid_end-2015+1)*12);
ssp126.mid.swf     = ssp126.mid.flooded + ssp126.mid.fh2osfc;
ssp126.mid.fsat    = ssp126.fsat(:,(mid_start-2015)*12+1:(mid_end-2015+1)*12);
ssp126.mid.perch   = ssp126.perch(:,(mid_start-2015)*12+1:(mid_end-2015+1)*12);
ssp126.mid.zwt     = ssp126.zwt(:,(mid_start-2015)*12+1:(mid_end-2015+1)*12);
ssp126.mid.tsa     = ssp126.tsa(:,(mid_start-2015)*12+1:(mid_end-2015+1)*12)- 273.15;
ssp126.mid.rain    = ssp126.rain(:,(mid_start-2015)*12+1:(mid_end-2015+1)*12);
flux.ssp126.mid.fgr    = flux.ssp126.fgr(:,(mid_start-2015)*12+1:(mid_end-2015+1)*12);
flux.ssp126.mid.surftsoi = flux.ssp126.surftsoi(:,(mid_start-2015)*12+1:(mid_end-2015+1)*12);
flux.ssp126.mid.soiltemp = flux.ssp126.soiltemp(:,(mid_start-2015)*12+1:(mid_end-2015+1)*12);
ice.ssp126.mid.soilice = ice.ssp126.soilice(:,(mid_start-2015)*12+1:(mid_end-2015+1)*12);
ice.ssp126.mid.surfice = ice.ssp126.surfice(:,(mid_start-2015)*12+1:(mid_end-2015+1)*12);
ice.ssp126.mid.snow   = ice.ssp126.snow(:,(mid_start-2015)*12+1:(mid_end-2015+1)*12);
ice.ssp126.mid.fsno   = ice.ssp126.fsno(:,(mid_start-2015)*12+1:(mid_end-2015+1)*12);
ssp126.end.flooded = ssp126.flooded(:,(end_start-2015)*12+1:(end_end-2015+1)*12);
ssp126.end.fh2osfc = ssp126.fh2osfc(:,(end_start-2015)*12+1:(end_end-2015+1)*12);
ssp126.end.swf     = ssp126.end.flooded + ssp126.end.fh2osfc;
ssp126.end.fsat    = ssp126.fsat(:,(end_start-2015)*12+1:(end_end-2015+1)*12);
ssp126.end.perch   = ssp126.perch(:,(end_start-2015)*12+1:(end_end-2015+1)*12);
ssp126.end.zwt     = ssp126.zwt(:,(end_start-2015)*12+1:(end_end-2015+1)*12);
ssp126.end.tsa     = ssp126.tsa(:,(end_start-2015)*12+1:(end_end-2015+1)*12)- 273.15;
ssp126.end.rain    = ssp126.rain(:,(end_start-2015)*12+1:(end_end-2015+1)*12);
flux.ssp126.end.fgr    = flux.ssp126.fgr(:,(end_start-2015)*12+1:(end_end-2015+1)*12);
flux.ssp126.end.surftsoi = flux.ssp126.surftsoi(:,(end_start-2015)*12+1:(end_end-2015+1)*12);
flux.ssp126.end.soiltemp = flux.ssp126.soiltemp(:,(end_start-2015)*12+1:(end_end-2015+1)*12);
ice.ssp126.end.soilice = ice.ssp126.soilice(:,(end_start-2015)*12+1:(end_end-2015+1)*12);
ice.ssp126.end.surfice = ice.ssp126.surfice(:,(end_start-2015)*12+1:(end_end-2015+1)*12);
ice.ssp126.end.snow   = ice.ssp126.snow(:,(end_start-2015)*12+1:(end_end-2015+1)*12);
ice.ssp126.end.fsno   = ice.ssp126.fsno(:,(end_start-2015)*12+1:(end_end-2015+1)*12);

ssp585.mid.flooded = ssp585.flooded(:,(mid_start-2015)*12+1:(mid_end-2015+1)*12);
ssp585.mid.fh2osfc = ssp585.fh2osfc(:,(mid_start-2015)*12+1:(mid_end-2015+1)*12);
ssp585.mid.swf     = ssp585.mid.flooded + ssp585.mid.fh2osfc;
ssp585.mid.fsat    = ssp585.fsat(:,(mid_start-2015)*12+1:(mid_end-2015+1)*12);
ssp585.mid.perch   = ssp585.perch(:,(mid_start-2015)*12+1:(mid_end-2015+1)*12);
ssp585.mid.zwt     = ssp585.zwt(:,(mid_start-2015)*12+1:(mid_end-2015+1)*12);
ssp585.mid.tsa     = ssp585.tsa(:,(mid_start-2015)*12+1:(mid_end-2015+1)*12) - 273.15;
ssp585.mid.rain    = ssp585.rain(:,(mid_start-2015)*12+1:(mid_end-2015+1)*12);
flux.ssp585.mid.fgr    = flux.ssp585.fgr(:,(mid_start-2015)*12+1:(mid_end-2015+1)*12);
flux.ssp585.mid.surftsoi = flux.ssp585.surftsoi(:,(mid_start-2015)*12+1:(mid_end-2015+1)*12);
flux.ssp585.mid.soiltemp = flux.ssp585.soiltemp(:,(mid_start-2015)*12+1:(mid_end-2015+1)*12);
ice.ssp585.mid.soilice = ice.ssp585.soilice(:,(mid_start-2015)*12+1:(mid_end-2015+1)*12);
ice.ssp585.mid.surfice = ice.ssp585.surfice(:,(mid_start-2015)*12+1:(mid_end-2015+1)*12);
ice.ssp585.mid.snow   = ice.ssp585.snow(:,(mid_start-2015)*12+1:(mid_end-2015+1)*12);
ice.ssp585.mid.fsno   = ice.ssp585.fsno(:,(mid_start-2015)*12+1:(mid_end-2015+1)*12);
ssp585.end.flooded = ssp585.flooded(:,(end_start-2015)*12+1:(end_end-2015+1)*12);
ssp585.end.fh2osfc = ssp585.fh2osfc(:,(end_start-2015)*12+1:(end_end-2015+1)*12);
ssp585.end.swf     = ssp585.end.flooded + ssp585.end.fh2osfc;
ssp585.end.fsat    = ssp585.fsat(:,(end_start-2015)*12+1:(end_end-2015+1)*12);
ssp585.end.perch   = ssp585.perch(:,(end_start-2015)*12+1:(end_end-2015+1)*12);
ssp585.end.zwt     = ssp585.zwt(:,(end_start-2015)*12+1:(end_end-2015+1)*12);
ssp585.end.tsa     = ssp585.tsa(:,(end_start-2015)*12+1:(end_end-2015+1)*12)- 273.15;
ssp585.end.rain    = ssp585.rain(:,(end_start-2015)*12+1:(end_end-2015+1)*12);
flux.ssp585.end.fgr    = flux.ssp585.fgr(:,(end_start-2015)*12+1:(end_end-2015+1)*12);
flux.ssp585.end.surftsoi = flux.ssp585.surftsoi(:,(end_start-2015)*12+1:(end_end-2015+1)*12);
flux.ssp585.end.soiltemp = flux.ssp585.soiltemp(:,(end_start-2015)*12+1:(end_end-2015+1)*12);
ice.ssp585.end.soilice = ice.ssp585.soilice(:,(end_start-2015)*12+1:(end_end-2015+1)*12);
ice.ssp585.end.surfice = ice.ssp585.surfice(:,(end_start-2015)*12+1:(end_end-2015+1)*12);
ice.ssp585.end.snow   = ice.ssp585.snow(:,(end_start-2015)*12+1:(end_end-2015+1)*12);
ice.ssp585.end.fsno   = ice.ssp585.fsno(:,(end_start-2015)*12+1:(end_end-2015+1)*12);

ssp126.mid.swf(lakein,:)   = NaN; ssp126.end.swf(lakein,:)   = NaN;
ssp126.mid.fsat(lakein,:)  = NaN; ssp126.end.fsat(lakein,:)  = NaN;
ssp126.mid.perch(lakein,:) = NaN; ssp126.end.perch(lakein,:) = NaN;
ssp126.mid.zwt(lakein,:)   = NaN; ssp126.end.zwt(lakein,:)   = NaN;
ssp126.mid.tsa(lakein,:)   = NaN; ssp126.end.tsa(lakein,:)   = NaN;
ssp126.mid.rain(lakein,:)  = NaN; ssp126.end.rain(lakein,:)  = NaN;
ssp585.mid.swf(lakein,:)   = NaN; ssp585.end.swf(lakein,:)   = NaN;
ssp585.mid.fsat(lakein,:)  = NaN; ssp585.end.fsat(lakein,:)  = NaN;
ssp585.mid.perch(lakein,:) = NaN; ssp585.end.perch(lakein,:) = NaN;
ssp585.mid.zwt(lakein,:)   = NaN; ssp585.end.zwt(lakein,:)   = NaN;
ssp585.mid.tsa(lakein,:)   = NaN; ssp585.end.tsa(lakein,:)   = NaN;
ssp585.mid.rain(lakein,:)  = NaN; ssp585.end.rain(lakein,:)  = NaN;
flux.ssp126.mid.fgr(lakein,:)    = NaN; flux.ssp126.end.fgr(lakein,:)    = NaN;
flux.ssp585.mid.fgr(lakein,:)    = NaN; flux.ssp585.end.fgr(lakein,:)    = NaN;
flux.ssp126.mid.surftsoi(lakein,:) = NaN; flux.ssp126.end.surftsoi(lakein,:) = NaN;
flux.ssp585.mid.surftsoi(lakein,:) = NaN; flux.ssp585.end.surftsoi(lakein,:) = NaN;
flux.ssp126.mid.soiltemp(lakein,:) = NaN; flux.ssp126.end.soiltemp(lakein,:) = NaN;
flux.ssp585.mid.soiltemp(lakein,:) = NaN; flux.ssp585.end.soiltemp(lakein,:) = NaN;
ice.ssp126.mid.soilice(lakein,:) = NaN; ice.ssp126.end.soilice(lakein,:) = NaN;
ice.ssp585.mid.soilice(lakein,:) = NaN; ice.ssp585.end.soilice(lakein,:) = NaN;
ice.ssp126.mid.surfice(lakein,:) = NaN; ice.ssp126.end.surfice(lakein,:) = NaN;
ice.ssp585.mid.surfice(lakein,:) = NaN; ice.ssp585.end.surfice(lakein,:) = NaN;
ice.ssp126.mid.snow(lakein,:) = NaN; ice.ssp126.end.snow(lakein,:) = NaN;
ice.ssp585.mid.snow(lakein,:) = NaN; ice.ssp585.end.snow(lakein,:) = NaN;
ice.ssp126.mid.fsno(lakein,:) = NaN; ice.ssp126.end.fsno(lakein,:) = NaN;
ice.ssp585.mid.fsno(lakein,:) = NaN; ice.ssp585.end.fsno(lakein,:) = NaN;

load(['../data/fsat_cal_12.mat'],'fsat_sea_cal');
irm = find(isnan(nanmean(fsat_sea_cal,2)));
historical.fsat(irm,:) = 0;
ssp126.mid.fsat(irm,:) = 0;
ssp126.end.fsat(irm,:) = 0;
ssp585.mid.fsat(irm,:) = 0;
ssp585.end.fsat(irm,:) = 0;

cmap = flipud(blue2red(17));

load('coastlines.mat');

% [fig0,axs0,cb0,pos0] = plot_two_axes(xv,yv, ...
%                     -nanmean(ssp126.end.zwt,2) + nanmean(historical.zwt,2),...
%                     -nanmean(ssp585.end.zwt,2) + nanmean(historical.zwt,2),...
%                     -2,2,'[m]');
% 
% [fig1,axs1,cb1,pos1] = plot_two_axes(xv,yv, ...
%                     -nanmean(ssp126.end.perch,2) + nanmean(historical.perch,2),...
%                     -nanmean(ssp585.end.perch,2) + nanmean(historical.perch,2),...
%                     -1,1,'[m]');
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

[fig2,axs2,cb2,pos2] = plot_two_axes(xv,yv, ...
                    (nanmean(ssp126.end.swf,2) - nanmean(historical.swf,2)).*area,...
                    (nanmean(ssp585.end.swf,2) - nanmean(historical.swf,2)).*area,...
                    -30.5, 30.5,'[km^2]');
colormap(axs2(1),cmap_new);
colormap(axs2(2),cmap_new);
add_title(axs2(1),'SSP126',20,'out');
add_title(axs2(2),'SSP585',20,'out');

ctl2d        = NaN(720,360);
ctl2d(index_lnd)        =  nanmean(historical.swf,2).*area;
foc2d_ssp126 = NaN(720,360);
foc2d_ssp126(index_lnd) = (nanmean(ssp126.end.swf,2) - nanmean(historical.swf,2)).*area; % [km^2]
foc2d_ssp585 = NaN(720,360);
foc2d_ssp585(index_lnd) = (nanmean(ssp585.end.swf,2) - nanmean(historical.swf,2)).*area; % [km^2]

foc1d_ssp126 = (nanmean(ssp126.end.swf,2) - nanmean(historical.swf,2)).*area; % [km^2]
foc1d_ssp585 = (nanmean(ssp585.end.swf,2) - nanmean(historical.swf,2)).*area; % [km^2]
axs2(3) = axes('Position',[pos2(1) + pos2(3) + 0.02 pos2(2) 0.25 axs2(1).Position(4)+axs2(1).Position(2)-axs2(2).Position(2)]);

ctl        = NaN(720,360);
ssp126_end = NaN(720,360);
ssp585_end = NaN(720,360);
ctl(index_lnd) = nanmean(historical.swf,2) .* area;
ssp126_end(index_lnd) = nanmean(ssp126.end.swf,2) .* area;

ctl2d(isnan(gladannual))        = NaN;
foc2d_ssp126(isnan(gladannual)) = NaN;
foc2d_ssp585(isnan(gladannual)) = NaN;
plot(axs2(3),movmean(nansum(foc2d_ssp126,1),5),-89.75:0.5:89.75,'b-', 'LineWidth',2);hold on; grid on;
plot(axs2(3),movmean(nansum(foc2d_ssp585,1),5),-89.75:0.5:89.75,'-','LineWidth',2); 
plot(axs2(3),zeros(360,1),-89.75:0.5:89.75,'k--','LineWidth',1); 
ylim([-60 80]);

set(axs2(1),'Color',[0.8 0.8 0.8]);
set(axs2(2),'Color',[0.8 0.8 0.8]);

legend('SSP126','SSP585');

[fig3,axs3,cb3,pos3] = plot_two_axes(xv,yv, ...
                    (nanmean(ssp126.end.fsat,2) - nanmean(historical.fsat,2)).*area,...
                    (nanmean(ssp585.end.fsat,2) - nanmean(historical.fsat,2)).*area,...
                    -120.5, 120.5,'[km^2]');
colormap(axs3(1),cmap_new);
colormap(axs3(2),cmap_new);
add_title(axs3(1),'SSP126',20,'out');
add_title(axs3(2),'SSP585',20,'out');

ctl2d        = NaN(720,360);
ctl2d(index_lnd)        =  nanmean(historical.fsat,2).*area;
foc2d_ssp126 = NaN(720,360);
foc2d_ssp126(index_lnd) = (nanmean(ssp126.end.fsat,2) - nanmean(historical.fsat,2)).*area; % [km^2]
foc2d_ssp585 = NaN(720,360);
foc2d_ssp585(index_lnd) = (nanmean(ssp585.end.fsat,2) - nanmean(historical.fsat,2)).*area; % [km^2]

foc1d_ssp126 = (nanmean(ssp126.end.fsat,2) - nanmean(historical.fsat,2)).*area; % [km^2]
foc1d_ssp585 = (nanmean(ssp585.end.fsat,2) - nanmean(historical.fsat,2)).*area; % [km^2]
axs3(3) = axes('Position',[pos2(1) + pos2(3) + 0.02 pos2(2) 0.25 axs3(1).Position(4)+axs3(1).Position(2)-axs3(2).Position(2)]);

% ctl        = NaN(720,360);
% ssp126_end = NaN(720,360);
% ssp585_end = NaN(720,360);
% ctl(index_lnd) = nanmean(historical.fsat,2) .* area;
% ssp126_end(index_lnd) = nanmean(ssp126.end.fsat,2) .* area;

ctl2d(isnan(gladannual))        = NaN;
foc2d_ssp126(isnan(gladannual)) = NaN;
foc2d_ssp585(isnan(gladannual)) = NaN;
plot(axs3(3),movmean(nansum(foc2d_ssp126,1),5)./movmean(nansum(ctl2d,1),5),-89.75:0.5:89.75,'b-', 'LineWidth',2);hold on; grid on;
plot(axs3(3),movmean(nansum(foc2d_ssp585,1),5)./movmean(nansum(ctl2d,1),5),-89.75:0.5:89.75,'-','LineWidth',2); 
plot(axs3(3),zeros(360,1),-89.75:0.5:89.75,'k--','LineWidth',1); 
ylim([-60 80]);

set(axs3(1),'Color',[0.8 0.8 0.8]);
set(axs3(2),'Color',[0.8 0.8 0.8]);

legend('SSP126','SSP585');

% [fig3,axs3,cb3,pos3] = plot_two_axes(xv,yv,    ...
%                     nanmean(ssp126.end.tsa,2) - nanmean(historical.tsa,2), ...
%                     nanmean(ssp585.end.tsa,2) - nanmean(historical.tsa,2), ...
%                     0,10,'[^{\circ}C]');
% colormap(jet(10));

% [fig4,axs4,cb4,pos4] = plot_two_axes(xv,yv, ...
%                     nanmean(flux.ssp126.end.soiltemp(:,ind),2)-273.15 ,...
%                     nanmean(flux.ssp585.end.soiltemp(:,ind),2)-273.15 ,...
%                     -10, 10,'[km^2]');
% colormap(axs4(1),cmap_new);
% colormap(axs4(2),cmap_new);
% add_title(axs4(1),'SSP126',20,'out');
% add_title(axs4(2),'SSP585',20,'out');

figure;

ind = [ [12:12:360] [1:12:360] [2:12:360] [3:12:360] ];

a = nanmean(flux.historical.soiltemp(:,ind),2);
b = nanmean(flux.ssp126.end.soiltemp(:,ind),2);
c = nanmean(flux.ssp585.end.soiltemp(:,ind),2);
d = (nanmean(ssp126.end.swf,2) - nanmean(historical.swf,2)).*area;
e = (nanmean(ssp585.end.swf,2) - nanmean(historical.swf,2)).*area;

aa = nanmean(ice.historical.soilice,2);
bb = nanmean(ice.ssp126.end.soilice,2);
cc = nanmean(ice.ssp585.end.soilice,2);
aa2 = nanmean(ice.historical.surfice,2);
bb2 = nanmean(ice.ssp126.end.surfice,2);
cc2 = nanmean(ice.ssp585.end.surfice,2);

aaa = nanmean(historical.rain ,2).*86400;
bbb = nanmean(ssp126.end.rain ,2).*86400;
ccc = nanmean(ssp585.end.rain ,2).*86400;

aaaa = nanmean(historical.perch,2);
bbbb = nanmean(ssp126.end.perch,2);
cccc = nanmean(ssp585.end.perch,2);

d(yc <= 50) = NaN; e(yc <= 50) = NaN;
ind3 = find(d > 0); ind4 = find(d < -0);
ind1 = find(e > 0); ind2 = find(e < -0);

% tmp1 = nanmean(reshape(nanmean(ssp585.end.swf(ind1,:),1)',[12 30]),2) - nanmean(reshape(nanmean(historical.swf(ind1,:),1)',[12 30]),2);
% tmp2 = nanmean(reshape(nanmean(ssp585.end.swf(ind2,:),1)',[12 30]),2) - nanmean(reshape(nanmean(historical.swf(ind2,:),1)',[12 30]),2);

figure; set(gcf,'Position',[10 10 800 400])
subplot(1,2,1);
tmp1 = nanmean(reshape(nanmean(historical.rain(ind1,:),1)',[12 30]),2);
tmp2 = nanmean(reshape(nanmean(historical.rain(ind2,:),1)',[12 30]),2);
plot(tmp1.*86400.*30,'b-','LineWidth',2); hold on; grid on;
plot(tmp2.*86400.*30,'r-','LineWidth',2);
set(gca,'FontSize',13);
ylabel('Rain [mm/month]','FontSize',15,'FontWeight','bold');
legend('Incrasing surface water regions','Decreasing surface water regions')
xlim([1 12]); ylim([0 80])
title('Control period','FontSize',20,'FontWeight','bold');

subplot(1,2,2);
tmp1 = nanmean(reshape(nanmean(ssp585.end.tsa(ind1,:),1)',[12 30]),2);
tmp2 = nanmean(reshape(nanmean(ssp585.end.tsa(ind2,:),1)',[12 30]),2);
plot(tmp1.*86400.*30,'b-','LineWidth',2); hold on; grid on;
plot(tmp2.*86400.*30,'r-','LineWidth',2);
set(gca,'FontSize',13);
ylabel('Rain [mm/month]','FontSize',15,'FontWeight','bold');
xlim([1 12]); 
title('Future period (SSP585)','FontSize',20,'FontWeight','bold');

figure; set(gcf,'Position',[10 10 800 400])
subplot(1,2,1);
tmp1 = nanmean(reshape(nanmean(historical.tsa(ind1,:),1)',[12 30]),2);
tmp2 = nanmean(reshape(nanmean(historical.tsa(ind2,:),1)',[12 30]),2);
plot(tmp1,'b-','LineWidth',2); hold on; grid on;
plot(tmp2,'r-','LineWidth',2);
set(gca,'FontSize',13);
ylabel('Rain [mm/month]','FontSize',15,'FontWeight','bold');
legend('Incrasing surface water regions','Decreasing surface water regions')
xlim([1 12]); 
title('Control period','FontSize',20,'FontWeight','bold');

subplot(1,2,2);
tmp1 = nanmean(reshape(nanmean(ssp585.end.tsa(ind1,:),1)',[12 30]),2);
tmp2 = nanmean(reshape(nanmean(ssp585.end.tsa(ind2,:),1)',[12 30]),2);
plot(tmp1,'b-','LineWidth',2); hold on; grid on;
plot(tmp2,'r-','LineWidth',2);
set(gca,'FontSize',13);
ylabel('Rain [mm/month]','FontSize',15,'FontWeight','bold');
xlim([1 12]); 
title('Future period (SSP585)','FontSize',20,'FontWeight','bold');


tmp1 = nanmean(reshape(nanmean(ssp585.end.tsa(ind1,:),1)',[12 30]),2) - nanmean(reshape(nanmean(historical.tsa(ind1,:),1)',[12 30]),2);
tmp2 = nanmean(reshape(nanmean(ssp585.end.tsa(ind2,:),1)',[12 30]),2) - nanmean(reshape(nanmean(historical.tsa(ind2,:),1)',[12 30]),2);
figure; set(gcf,'Position',[10 10 1200 1000])
subplot(2,2,1);
plot(tmp1,'b-','LineWidth',2); hold on; grid on;
plot(tmp2,'r-','LineWidth',2);
set(gca,'FontSize',13);
ylabel('Temperature Change [C]','FontSize',15,'FontWeight','bold');
xlim([1 12]); 
legend('Incrasing surface water regions','Decreasing surface water regions')

tmp1 = nanmean(reshape(nanmean(ssp585.end.rain(ind1,:),1)',[12 30]),2) - nanmean(reshape(nanmean(historical.rain(ind1,:),1)',[12 30]),2);
tmp2 = nanmean(reshape(nanmean(ssp585.end.rain(ind2,:),1)',[12 30]),2) - nanmean(reshape(nanmean(historical.rain(ind2,:),1)',[12 30]),2);
subplot(2,2,2);
plot(tmp1.*86400.*30,'b-','LineWidth',2); hold on; grid on;
plot(tmp2.*86400.*30,'r-','LineWidth',2);
set(gca,'FontSize',13);
ylabel('Rain Change [mm/month]','FontSize',15,'FontWeight','bold');
xlim([1 12])
% 
% tmp1 = nanmean(reshape(nanmean(ice.ssp126.end.soilice(ind1,:),1)',[12 30]),2) - nanmean(reshape(nanmean(ice.historical.soilice(ind1,:),1)',[12 30]),2);
% tmp2 = nanmean(reshape(nanmean(ice.ssp126.end.soilice(ind2,:),1)',[12 30]),2) - nanmean(reshape(nanmean(ice.historical.soilice(ind2,:),1)',[12 30]),2);
% 
tmp1 = nanmean(reshape(nanmean(ice.ssp585.end.soilice(ind1,:),1)',[12 30]),2) - nanmean(reshape(nanmean(ice.historical.soilice(ind1,:),1)',[12 30]),2);
tmp2 = nanmean(reshape(nanmean(ice.ssp585.end.soilice(ind2,:),1)',[12 30]),2) - nanmean(reshape(nanmean(ice.historical.soilice(ind2,:),1)',[12 30]),2);
subplot(2,2,3);
plot(tmp1,'b-','LineWidth',2); hold on; grid on;
plot(tmp2,'r-','LineWidth',2);
set(gca,'FontSize',13);
ylabel('Soil ice Change [kg/m^2]','FontSize',15,'FontWeight','bold');
xlim([1 12])
tmp1 = nanmean(reshape(nanmean(ice.ssp585.end.surfice(ind1,:),1)',[12 30]),2) - nanmean(reshape(nanmean(ice.historical.surfice(ind1,:),1)',[12 30]),2);
tmp2 = nanmean(reshape(nanmean(ice.ssp585.end.surfice(ind2,:),1)',[12 30]),2) - nanmean(reshape(nanmean(ice.historical.surfice(ind2,:),1)',[12 30]),2);
subplot(2,2,4);
plot(tmp1,'b-','LineWidth',2); hold on; grid on;
plot(tmp2,'r-','LineWidth',2);
set(gca,'FontSize',13);
ylabel('Surface ice Change [kg/m^2]','FontSize',15,'FontWeight','bold');
xlim([1 12])

figure;
g1 = [ones(length(ind1),1); ones(length(ind1),1).*2; ones(length(ind1),1).*3];
g2 = [ones(length(ind2),1).*4; ones(length(ind2),1).*5; ones(length(ind2),1).*6];
h(1) = boxchart(g1,[a(ind1); b(ind1); c(ind1)]-273.15,'MarkerStyle','none'); hold on; grid on;
h(2) = boxchart(g2,[a(ind2); b(ind2); c(ind2)]-273.15,'MarkerStyle','none');
set(h,{'linew'},{2});
set(gca,'FontSize',13);
set(gca,'XTick',[1 2 3 4 5 6]);
set(gca,'XTickLabel',{'CTL','SSP126','SSP585','CTL','SSP126','SSP585'});
legend({'Increasing region','Decreasing region'});
ylabel('Soil temperature [^{\circ}C]','FontSize',15,'FontWeight','bold');
ylim([-27 15]);

figure;
g1 = [ones(length(ind1),1); ones(length(ind1),1).*2; ones(length(ind1),1).*3];
g2 = [ones(length(ind2),1).*4; ones(length(ind2),1).*5; ones(length(ind2),1).*6];
h(1) = boxchart(g1,[aa(ind1); bb(ind1); cc(ind1)],'MarkerStyle','none'); hold on; grid on;
h(2) = boxchart(g2,[aa(ind2); bb(ind2); cc(ind2)],'MarkerStyle','none');
set(h,{'linew'},{2});
set(gca,'FontSize',13);
set(gca,'XTick',[1 2 3 4 5 6]);
set(gca,'XTickLabel',{'CTL','SSP126','SSP585','CTL','SSP126','SSP585'});
legend({'Increasing region','Decreasing region'});
ylabel('Soil ice [kg]','FontSize',15,'FontWeight','bold');

figure;
g1 = [ ones(length(ind1),1).*1; ones(length(ind1),1).*2];
g2 = [ ones(length(ind2),1).*3; ones(length(ind2),1).*4];

h(1) = boxchart(g1,[bbb(ind1) - aaa(ind1) ; ccc(ind1) - aaa(ind1)],'MarkerStyle','none'); hold on; grid on;
h(2) = boxchart(g2,[bbb(ind2) - aaa(ind2) ; ccc(ind2) - aaa(ind2)],'MarkerStyle','none');

set(h,{'linew'},{2});
set(gca,'FontSize',13);
set(gca,'XTick',[1 2 3 4]);
set(gca,'XTickLabel',{'SSP126','SSP585','SSP126','SSP585'});
legend({'Increasing region','Decreasing region'});
ylabel('Change of rain + snowmelt [mm/year]','FontSize',15,'FontWeight','bold');
ylim([-0.2 1]);

figure;
g1 = [ ones(length(ind1),1).*1; ones(length(ind1),1).*2];
g2 = [ ones(length(ind2),1).*3; ones(length(ind2),1).*4];

h(1) = boxchart(g1,[b(ind1) - a(ind1); c(ind1) - a(ind1)],'MarkerStyle','none'); hold on; grid on;
h(2) = boxchart(g2,[b(ind2) - a(ind2); c(ind2) - a(ind2)],'MarkerStyle','none');

set(h,{'linew'},{2});
set(gca,'FontSize',13);
set(gca,'XTick',[1 2 3 4]);
set(gca,'XTickLabel',{'SSP126','SSP585','SSP126','SSP585'});
legend({'Increasing region','Decreasing region'});
ylabel('Change of surface temperature [kg]','FontSize',15,'FontWeight','bold');

figure;
g1 = [ ones(length(ind1),1).*1; ones(length(ind1),1).*2];
g2 = [ ones(length(ind2),1).*3; ones(length(ind2),1).*4];

h(1) = boxchart(g1,[bb(ind1) - aa(ind1); cc(ind1) - aa(ind1)],'MarkerStyle','none'); hold on; grid on;
h(2) = boxchart(g2,[bb(ind2) - aa(ind2); cc(ind2) - aa(ind2)],'MarkerStyle','none');

set(h,{'linew'},{2});
set(gca,'FontSize',13);
set(gca,'XTick',[1 2 3 4]);
set(gca,'XTickLabel',{'SSP126','SSP585','SSP126','SSP585'});
legend({'Increasing region','Decreasing region'});
ylabel('Change of soil ice [kg]','FontSize',15,'FontWeight','bold');

figure;
g1 = [ ones(length(ind1),1).*1; ones(length(ind1),1).*2];
g2 = [ ones(length(ind2),1).*3; ones(length(ind2),1).*4];

h(1) = boxchart(g1,[bb2(ind1) - aa2(ind1); cc2(ind1) - aa2(ind1)],'MarkerStyle','none'); hold on; grid on;
h(2) = boxchart(g2,[bb2(ind2) - aa2(ind2); cc2(ind2) - aa2(ind2)],'MarkerStyle','none');

set(h,{'linew'},{2});
set(gca,'FontSize',13);
set(gca,'XTick',[1 2 3 4]);
set(gca,'XTickLabel',{'SSP126','SSP585','SSP126','SSP585'});
legend({'Increasing region','Decreasing region'});
ylabel('Change of surface ice [kg]','FontSize',15,'FontWeight','bold');

figure;
g1 = [ ones(length(ind1),1).*1; ones(length(ind1),1).*2];
g2 = [ ones(length(ind2),1).*3; ones(length(ind2),1).*4];

h(1) = boxchart(g1,[bbbb(ind1) - aaaa(ind1); cccc(ind1) - aaaa(ind1)],'MarkerStyle','none'); hold on; grid on;
h(2) = boxchart(g2,[bbbb(ind2) - aaaa(ind2); cccc(ind2) - aaaa(ind2)],'MarkerStyle','none');

set(h,{'linew'},{2});
set(gca,'FontSize',13);
set(gca,'XTick',[1 2 3 4]);
set(gca,'XTickLabel',{'SSP126','SSP585','SSP126','SSP585'});
legend({'Increasing region','Decreasing region'});
ylabel('Change of perched water table [m]','FontSize',15,'FontWeight','bold');
ylim([-0.4 1.2]);

[fig4,axs4,cb4,pos4] = plot_two_axes(xv,yv,    ...
                    nanmean(ice.ssp126.end.soilice,2) - nanmean(ice.historical.soilice,2), ...
                    nanmean(ice.ssp585.end.soilice,2) - nanmean(ice.historical.soilice,2), ...
                    -30,30,'[-]');

% plot([0.5 4.5],[0 0],'r--','LineWidth',2);

figure;
histogram(a(ind)); hold on;
histogram(b(ind));
histogram(c(ind));
% figure;
% plot(nanmean(ssp585.end.perch,2) - nanmean(historical.perch,2),(nanmean(ssp585.end.swf,2) - nanmean(historical.swf,2)).*area,'bx');
PCT_NAT_PFT = ncread('../inputdata/landuse.timeseries_0.5x0.5_SSP5_RCP85_simyr2015-2100_GLOBAL_1d_c240209.nc','PCT_NAT_PFT');
PCT_NAT_PFT = nanmean(PCT_NAT_PFT,3);
MONTHLY_LAI = ncread('../inputdata/surfdata_GLOBE_1d_calibrated.nc','MONTHLY_LAI');
MONTHLY_SAI = ncread('../inputdata/surfdata_GLOBE_1d_calibrated.nc','MONTHLY_SAI');
STD_ELEV    = ncread('../inputdata/surfdata_GLOBE_1d_calibrated.nc','STD_ELEV');
TOPO        = ncread('../inputdata/surfdata_GLOBE_1d_calibrated.nc','TOPO');

lai = nansum(repmat(PCT_NAT_PFT,1,1,12) .* MONTHLY_LAI,2) ./100;
lai = reshape(lai,[66924,12]);
sai = nansum(repmat(PCT_NAT_PFT,1,1,12) .* MONTHLY_SAI,2) ./100;
sai = reshape(sai,[66924,12]);

ind = find(yc >= 40);
pft = {'not_vegetated                      '; ...
       'needleleaf_evergreen_temperate_tree'; ...
       'needleleaf_evergreen_boreal_tree   '; ...
       'needleleaf_deciduous_boreal_tree   '; ...
       'broadleaf_evergreen_tropical_tree  '; ...
       'broadleaf_evergreen_temperate_tree '; ...
       'broadleaf_deciduous_tropical_tree  '; ...
       'broadleaf_deciduous_temperate_tree '; ...
       'broadleaf_deciduous_boreal_tree    '; ...
       'broadleaf_evergreen_shrub          '; ...
       'broadleaf_deciduous_temperate_shrub'; ...
       'broadleaf_deciduous_boreal_shrub   '; ...
       'c3_arctic_grass                    '; ...
       'c3_non-arctic_grass                '; ...
       'c4_grass                           '; ...
       'c3_crop                            '; ...
       'c3_irrigated                       '};

figure;
k = 1;
for i = [1 3 4 12 13]
    subplot(5,1,k);
    patch(xv(:,ind),yv(:,ind),PCT_NAT_PFT(ind,i),'LineStyle','none'); colorbar;
    title(pft{i},'Interpreter','none')
    k = k + 1;
end
figure;
k = 1;
for i = [1 3 4 12 13]
    subplot(5,1,k);
    [f1,x1] = ksdensity(PCT_NAT_PFT(ind1,i),'Function','cdf');
    [f2,x2] = ksdensity(PCT_NAT_PFT(ind2,i),'Function','cdf');
    plot(x1,f1,'b-','LineWidth',2); hold on; grid on;
    plot(x2,f2,'r-','LineWidth',2);
    title(pft{i},'Interpreter','none');
    k = k + 1;
end


figure;
patch(xv,yv,nanmean(lai,2),'LineStyle','none'); hold on;

figure;
plot(nanmean(PCT_NAT_PFT(ind1,:),1),'b-','LineWidth',2); hold on; grid on;
plot(nanmean(PCT_NAT_PFT(ind2,:),1),'r-','LineWidth',2);

ind = find(nanmean(PCT_NAT_PFT(:,16,:),3) > 10);
ind1 = [ [5:12:360] [6:12:360] [7:12:360] [8:12:360] ];
y1 = nanmean(historical.swf(:,ind1),2).*area;
y2 = nanmean(ssp585.end.swf(:,ind1),2).*area;
x  = nanmean(PCT_NAT_PFT(:,16,:),3)./100.*area;

y1 = y1(ind);
y2 = y2(ind);
x  = x(ind);
f1 = y1./x; f1(f1 > 1) = 1;
f2 = y2./x; f2(f2 > 1) = 1;

figure;
histogram(f2 - f1);

ind1 = find(e > 5); ind2 = find(e < -5);
b1 = nanmean(ssp585.end.rain(ind1,:),1) ;
b2 = nanmean(ssp585.end.rain(ind2,:),1) ;
a1 = nanmean(historical.rain(ind1,:),1) ;
a2 = nanmean(historical.rain(ind2,:),1) ;

figure;
plot(nanmean(reshape(a1',[12 30]),2),'b--','LineWidth',2); hold on;
plot(nanmean(reshape(a2',[12 30]),2),'r--','LineWidth',2);
plot(nanmean(reshape(b1',[12 30]),2),'b-','LineWidth',2); 
plot(nanmean(reshape(b2',[12 30]),2),'r-','LineWidth',2);

title('Soil ice','FontSize',15,'FontWeight','bold');
xlim([1 12]);
legend('Increasing region - Historical','Decreasing region - Historical','Increasing region - SSP585','Decreasing region - SSP585','FontSize',15,'FontWeight','bold');

b1 = nanmean(ice.ssp585.end.soilice(ind1,:),1) ;
b2 = nanmean(ice.ssp585.end.soilice(ind2,:),1) ;
a1 = nanmean(ice.historical.soilice(ind1,:),1) ;
a2 = nanmean(ice.historical.soilice(ind2,:),1) ;

figure;
plot(nanmean(reshape(a1',[12 30]),2),'b--','LineWidth',2); hold on;
plot(nanmean(reshape(a2',[12 30]),2),'r--','LineWidth',2);
plot(nanmean(reshape(b1',[12 30]),2),'b-','LineWidth',2); 
plot(nanmean(reshape(b2',[12 30]),2),'r-','LineWidth',2);

title('Soil ice','FontSize',15,'FontWeight','bold');
xlim([1 12]);
legend('Increasing region - Historical','Decreasing region - Historical','Increasing region - SSP585','Decreasing region - SSP585','FontSize',15,'FontWeight','bold');


b1 = nanmean(flux.ssp585.end.soiltemp(ind1,:)-273.15,1) ;
b2 = nanmean(flux.ssp585.end.soiltemp(ind2,:)-273.15,1) ;
a1 = nanmean(flux.historical.soiltemp(ind1,:)-273.15,1) ;
a2 = nanmean(flux.historical.soiltemp(ind2,:)-273.15,1) ;
figure;
plot(nanmean(reshape(a1',[12 30]),2),'b--','LineWidth',2); hold on;
plot(nanmean(reshape(a2',[12 30]),2),'r--','LineWidth',2);
plot(nanmean(reshape(b1',[12 30]),2),'b-','LineWidth',2); 
plot(nanmean(reshape(b2',[12 30]),2),'r-','LineWidth',2);


title('Soil temperature','FontSize',15,'FontWeight','bold');
xlim([1 12]);
legend('Increasing region - Historical','Decreasing region - Historical','Increasing region - SSP585','Decreasing region - SSP585','FontSize',15,'FontWeight','bold');

b1 = nanmean(flux.ssp585.end.fgr(ind1,:),1) ;
b2 = nanmean(flux.ssp585.end.fgr(ind2,:),1) ;
a1 = nanmean(flux.historical.fgr(ind1,:),1) ;
a2 = nanmean(flux.historical.fgr(ind2,:),1) ;
figure;
plot(nanmean(reshape(a1',[12 30]),2),'b--','LineWidth',2); hold on;
plot(nanmean(reshape(a2',[12 30]),2),'r--','LineWidth',2);
plot(nanmean(reshape(b1',[12 30]),2),'b-','LineWidth',2); 
plot(nanmean(reshape(b2',[12 30]),2),'r-','LineWidth',2);
title('Ground heat flux','FontSize',15,'FontWeight','bold');
xlim([1 12]);
legend('Increasing region - Historical','Decreasing region - Historical','Increasing region - SSP585','Decreasing region - SSP585','FontSize',15,'FontWeight','bold');

b1 = nanmean(ice.ssp585.end.fsno(ind1,:),1) ;
b2 = nanmean(ice.ssp585.end.fsno(ind2,:),1) ;
a1 = nanmean(ice.historical.fsno(ind1,:),1) ;
a2 = nanmean(ice.historical.fsno(ind2,:),1) ;
figure;
plot(nanmean(reshape(b1',[12 30]),2),'b-','LineWidth',2); hold on;
plot(nanmean(reshape(b2',[12 30]),2),'r-','LineWidth',2);
plot(nanmean(reshape(a1',[12 30]),2),'b--','LineWidth',2); 
plot(nanmean(reshape(a2',[12 30]),2),'r--','LineWidth',2);
title('Snow fraction','FontSize',15,'FontWeight','bold');
xlim([1 12]);
legend('Increasing region - Historical','Decreasing region - Historical','Increasing region - SSP585','Decreasing region - SSP585','FontSize',15,'FontWeight','bold');

figure;
subplot(1,2,1);
cdfplot(TOPO(ind1)); hold on;
cdfplot(TOPO(ind2)); grid on;
title('Mean Elevation [m]','FontSize',15,'FontWeight','bold');
legend('Increasing region','Decreasing region','FontSize',15,'FontWeight','bold')
subplot(1,2,2);
cdfplot(STD_ELEV(ind1)); hold on;
cdfplot(STD_ELEV(ind2)); grid on;
title('Standard Deviation [m]','FontSize',15,'FontWeight','bold');