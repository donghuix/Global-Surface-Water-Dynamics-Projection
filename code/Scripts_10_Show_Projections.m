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
xv   = ncread('../data/domain_lnd_GLOBE_1d.nc','xv');
yv   = ncread('../data/domain_lnd_GLOBE_1d.nc','yv');
xc   = ncread('../data/domain_lnd_GLOBE_1d.nc','xc');
yc   = ncread('../data/domain_lnd_GLOBE_1d.nc','yc');
area = ncread('../data/domain_lnd_GLOBE_1d.nc','area');
[lakein,lakein2d] = getLakeIndex(e3sm_input);

historical      = load(['../projection/projection_BC_' model '_historical.mat'],     'rain','fh2osfc','flooded');
ssp126          = load(['../projection/projection_BC_' model '_ssp126.mat'],         'rain','fh2osfc','flooded');
ssp585          = load(['../projection/projection_BC_' model '_ssp585.mat'],         'rain','fh2osfc','flooded');
ice.historical  = load(['../projection/projection_BC_' model '_historical_ice.mat'], 'qsnow','snow');
ice.ssp126      = load(['../projection/projection_BC_' model '_ssp126_ice.mat'],     'qsnow','snow');
ice.ssp585      = load(['../projection/projection_BC_' model '_ssp585_ice.mat'],     'qsnow','snow');
flux.historical = load(['../projection/projection_BC_' model '_historical_flux.mat'],'latent');
flux.ssp126     = load(['../projection/projection_BC_' model '_ssp126_flux.mat'],    'latent');
flux.ssp585     = load(['../projection/projection_BC_' model '_ssp585_flux.mat'],    'latent');

historical.rain        = historical.rain(:,1:360);
historical.swf         = historical.fh2osfc(:,1:360)+historical.flooded(:,1:360);
ice.historical.qsnow   = ice.historical.qsnow(:,1:360);
ice.historical.snow    = ice.historical.snow(:,1:360);
historical.prec        = (historical.rain + ice.historical.snow).*86400;
historical.liquid      = (historical.rain + ice.historical.qsnow).*86400;
flux.historical.latent = flux.historical.latent(:,1:360);

ssp126.mid.rain        = ssp126.rain(:,(mid_start-2015)*12+1:(mid_end-2015+1)*12);
ssp126.mid.swf         = ssp126.fh2osfc(:,(mid_start-2015)*12+1:(mid_end-2015+1)*12) + ...
                         ssp126.flooded(:,(mid_start-2015)*12+1:(mid_end-2015+1)*12);
ice.ssp126.mid.qsnow   = ice.ssp126.qsnow(:,(mid_start-2015)*12+1:(mid_end-2015+1)*12);
ice.ssp126.mid.snow    = ice.ssp126.snow(:,(mid_start-2015)*12+1:(mid_end-2015+1)*12);
ssp126.mid.prec        = ssp126.mid.rain + ice.ssp126.mid.snow;
ssp126.mid.liquid      = ssp126.mid.rain + ice.ssp126.mid.qsnow;
flux.ssp126.mid.latent = flux.ssp126.latent(:,(mid_start-2015)*12+1:(mid_end-2015+1)*12);

ssp126.end.rain        = ssp126.rain(:,(end_start-2015)*12+1:(end_end-2015+1)*12);
ssp126.end.swf         = ssp126.fh2osfc(:,(end_start-2015)*12+1:(end_end-2015+1)*12) + ...
                         ssp126.flooded(:,(end_start-2015)*12+1:(end_end-2015+1)*12);
ice.ssp126.end.qsnow   = ice.ssp126.qsnow(:,(end_start-2015)*12+1:(end_end-2015+1)*12);
ice.ssp126.end.snow    = ice.ssp126.snow(:,(end_start-2015)*12+1:(end_end-2015+1)*12);
ssp126.end.prec        = (ssp126.end.rain + ice.ssp126.end.snow).*86400;
ssp126.end.liquid      = (ssp126.end.rain + ice.ssp126.end.qsnow).*86400;
flux.ssp126.end.latent = flux.ssp126.latent(:,(end_start-2015)*12+1:(end_end-2015+1)*12);

ssp585.mid.rain        = ssp585.rain(:,(mid_start-2015)*12+1:(mid_end-2015+1)*12);
ssp585.mid.swf         = ssp585.fh2osfc(:,(mid_start-2015)*12+1:(mid_end-2015+1)*12) + ...
                         ssp585.flooded(:,(mid_start-2015)*12+1:(mid_end-2015+1)*12);
ice.ssp585.mid.qsnow   = ice.ssp585.qsnow(:,(mid_start-2015)*12+1:(mid_end-2015+1)*12);
ice.ssp585.mid.snow    = ice.ssp585.snow(:,(mid_start-2015)*12+1:(mid_end-2015+1)*12);
ssp585.mid.prec        = ssp585.mid.rain + ice.ssp585.mid.snow;
ssp585.mid.liquid      = ssp585.mid.rain + ice.ssp585.mid.qsnow;
flux.ssp585.mid.latent = flux.ssp585.latent(:,(mid_start-2015)*12+1:(mid_end-2015+1)*12);

ssp585.end.rain        = ssp585.rain(:,(end_start-2015)*12+1:(end_end-2015+1)*12);
ssp585.end.swf         = ssp585.fh2osfc(:,(end_start-2015)*12+1:(end_end-2015+1)*12) + ...
                         ssp585.flooded(:,(end_start-2015)*12+1:(end_end-2015+1)*12);
ice.ssp585.end.qsnow   = ice.ssp585.qsnow(:,(end_start-2015)*12+1:(end_end-2015+1)*12);
ice.ssp585.end.snow    = ice.ssp585.snow(:,(end_start-2015)*12+1:(end_end-2015+1)*12);
ssp585.end.prec        = (ssp585.end.rain + ice.ssp585.end.snow).*86400;
ssp585.end.liquid      = (ssp585.end.rain + ice.ssp585.end.qsnow).*86400;
flux.ssp585.end.latent = flux.ssp585.latent(:,(end_start-2015)*12+1:(end_end-2015+1)*12);

% convert W/m^2 to mm/day
[ historical.ET ] = E1toE2( flux.historical.latent,2 );
[ ssp126.mid.ET ] = E1toE2( flux.ssp126.mid.latent,2 );
[ ssp126.end.ET ] = E1toE2( flux.ssp126.end.latent,2 );
[ ssp585.mid.ET ] = E1toE2( flux.ssp585.mid.latent,2 );
[ ssp585.end.ET ] = E1toE2( flux.ssp585.end.latent,2 );

cmap = flipud(blue2red(17));
historical.prec(lakein,:)   = NaN;
ssp126.mid.prec(lakein,:)   = NaN; ssp126.end.prec(lakein,:)   = NaN;
ssp585.mid.prec(lakein,:)   = NaN; ssp585.end.prec(lakein,:)   = NaN;
historical.liquid(lakein,:) = NaN;
ssp126.mid.liquid(lakein,:) = NaN; ssp126.end.liquid(lakein,:) = NaN;
ssp585.mid.liquid(lakein,:) = NaN; ssp585.end.liquid(lakein,:) = NaN;
historical.ET(lakein,:)     = NaN;
ssp126.mid.ET(lakein,:)     = NaN; ssp126.end.ET(lakein,:)     = NaN;
ssp585.mid.ET(lakein,:)     = NaN; ssp585.end.ET(lakein,:)     = NaN;
historical.swf(lakein,:)    = NaN;
ssp126.mid.swf(lakein,:)    = NaN; ssp126.end.swf(lakein,:)    = NaN;
ssp585.mid.swf(lakein,:)    = NaN; ssp585.end.swf(lakein,:)    = NaN;

historical.PE = historical.liquid - historical.ET;
ssp126.mid.PE = ssp126.mid.liquid - ssp126.mid.ET;
ssp126.end.PE = ssp126.end.liquid - ssp126.end.ET;
ssp585.mid.PE = ssp585.mid.liquid - ssp585.mid.ET;
ssp585.end.PE = ssp585.end.liquid - ssp585.end.ET;

cmap = getPanoply_cMap('UKM_hadcrut_10');
cmap(7,3) = 1;
cmap(12,:)= [];
cmap(6,:)= [];
cmap = flipud(cmap);

N = 1;
cmap_new = NaN(15*2*N+1,3);
cmap_new(1:5,:)   = repmat(cmap(1,:),5,1);
cmap_new(6:9,:)   = repmat(cmap(2,:),4,1);
cmap_new(10:12,:) = repmat(cmap(3,:),3,1);
cmap_new(13:14,:) = repmat(cmap(4,:),2,1);
cmap_new(15,:)    = repmat(cmap(5,:),1,1);
cmap_new(16,:)    = repmat(cmap(6,:),1,1);
cmap_new(17,:)    = repmat(cmap(7,:),1,1);
cmap_new(18:19,:) = repmat(cmap(8,:),2,1);
cmap_new(20:22,:) = repmat(cmap(9,:),3,1);
cmap_new(23:26,:) = repmat(cmap(10,:),4,1);
cmap_new(27:31,:) = repmat(cmap(11,:),5,1);

[fig1,axs1,cb1,pos1] = plot_two_axes(xv,yv, ...
                    nanmean(ssp126.end.liquid,2) - nanmean(historical.liquid,2),...
                    nanmean(ssp585.end.liquid,2) - nanmean(historical.liquid,2),...
                    -2,2,'[mm/day]'); colormap(cmap_new);
tks = [-2 -2+2/15.5*5 -2+2/15.5*9 -2+2/15.5*12 -2+2/15.5*14  0  2-2/15.5*14 2-2/15.5*12 2-2/15.5*9 2-2/15.5*5 2];
cb1.Ticks = tks;
cb1.TickLabels = string(round(tks,1));

% [fig2,axs2,cb2,pos2] = plot_two_axes(xv,yv, ...
%                     (nanmean(ssp126.end.prec,2) - nanmean(historical.prec,2))./nanmean(historical.prec,2).*100,...
%                     (nanmean(ssp585.end.prec,2) - nanmean(historical.prec,2))./nanmean(historical.prec,2).*100,...
%                     -31.875,31.875,'[%]');
% 
% [fig2,axs2,cb2,pos2] = plot_two_axes(xv,yv, ...
%                     (nanmean(ssp126.end.ET,2) - nanmean(historical.ET,2))./nanmean(historical.ET,2).*100,...
%                     (nanmean(ssp585.end.ET,2) - nanmean(historical.ET,2))./nanmean(historical.ET,2).*100,...
%                     -31.875,31.875,'[%]');

[fig3,axs3,cb3,pos3] = plot_two_axes(xv,yv, ...
                    nanmean(ssp126.end.PE,2) - nanmean(historical.PE,2),...
                    nanmean(ssp585.end.PE,2) - nanmean(historical.PE,2),...
                    -2,2,'[mm/day]'); colormap(cmap_new);
tks = [-2 -2+2/15.5*5 -2+2/15.5*9 -2+2/15.5*12 -2+2/15.5*14  0  2-2/15.5*14 2-2/15.5*12 2-2/15.5*9 2-2/15.5*5 2];
cb3.Ticks = tks;
cb3.TickLabels = string(round(tks,1));

[fig4,axs4,cb4,pos4] = plot_two_axes(xv,yv, ...
                    nanmean(ssp126.end.ET,2) - nanmean(historical.ET,2),...
                    nanmean(ssp585.end.ET,2) - nanmean(historical.ET,2),...
                    -2,2,'[mm/day]'); colormap(cmap_new);
tks = [-2 -2+2/15.5*5 -2+2/15.5*9 -2+2/15.5*12 -2+2/15.5*14  0  2-2/15.5*14 2-2/15.5*12 2-2/15.5*9 2-2/15.5*5 2];
cb4.Ticks = tks;
cb4.TickLabels = string(round(tks,1));

% x = nanmean(ssp585.end.PE,2)  - nanmean(historical.PE,2);
% y = nanmean(ssp585.end.swf,2) - nanmean(historical.swf,2);
% ind1 = find(x > 0 & y < 0);
% ind2 = find(x < 0 & y > 0);
% plot(axs3(2),xc(ind1),yc(ind1),'w.');
% plot(axs3(2),xc(ind2),yc(ind2),'g.');