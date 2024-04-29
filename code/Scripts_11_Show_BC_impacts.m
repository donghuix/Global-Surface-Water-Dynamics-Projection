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

historical = load(['../projection/projection_' model '_historical.mat'],'flooded','fh2osfc');
ssp126     = load(['../projection/projection_' model '_ssp126.mat'],'flooded','fh2osfc');
ssp585     = load(['../projection/projection_' model '_ssp585.mat'],'flooded','fh2osfc');

historical.bc = load(['../projection/projection_BC_' model '_historical.mat'],'flooded','fh2osfc');
ssp126.bc     = load(['../projection/projection_BC_' model '_ssp126.mat'],'flooded','fh2osfc');
ssp585.bc     = load(['../projection/projection_BC_' model '_ssp585.mat'],'flooded','fh2osfc');

historical.bc.flooded = historical.bc.flooded(:,1:360);
historical.bc.fh2osfc = historical.bc.fh2osfc(:,1:360);
historical.bc.swf     = historical.bc.flooded + historical.bc.fh2osfc;
ssp126.bc.mid.flooded = ssp126.bc.flooded(:,(mid_start-2015)*12+1:(mid_end-2015+1)*12);
ssp126.bc.mid.fh2osfc = ssp126.bc.fh2osfc(:,(mid_start-2015)*12+1:(mid_end-2015+1)*12);
ssp126.bc.mid.swf     = ssp126.bc.mid.flooded + ssp126.bc.mid.fh2osfc;
ssp126.bc.end.flooded = ssp126.bc.flooded(:,(end_start-2015)*12+1:(end_end-2015+1)*12);
ssp126.bc.end.fh2osfc = ssp126.bc.fh2osfc(:,(end_start-2015)*12+1:(end_end-2015+1)*12);
ssp126.bc.end.swf     = ssp126.bc.end.flooded + ssp126.bc.end.fh2osfc;
ssp585.bc.mid.flooded = ssp585.bc.flooded(:,(mid_start-2015)*12+1:(mid_end-2015+1)*12);
ssp585.bc.mid.fh2osfc = ssp585.bc.fh2osfc(:,(mid_start-2015)*12+1:(mid_end-2015+1)*12);
ssp585.bc.mid.swf     = ssp585.bc.mid.flooded + ssp585.bc.mid.fh2osfc;
ssp585.bc.end.flooded = ssp585.bc.flooded(:,(end_start-2015)*12+1:(end_end-2015+1)*12);
ssp585.bc.end.fh2osfc = ssp585.bc.fh2osfc(:,(end_start-2015)*12+1:(end_end-2015+1)*12);
ssp585.bc.end.swf     = ssp585.bc.end.flooded + ssp585.bc.end.fh2osfc;

historical.flooded = historical.flooded(:,1:360);
historical.fh2osfc = historical.fh2osfc(:,1:360);
historical.swf     = historical.flooded + historical.fh2osfc;
ssp126.mid.flooded = ssp126.flooded(:,(mid_start-2015)*12+1:(mid_end-2015+1)*12);
ssp126.mid.fh2osfc = ssp126.fh2osfc(:,(mid_start-2015)*12+1:(mid_end-2015+1)*12);
ssp126.mid.swf     = ssp126.mid.flooded + ssp126.mid.fh2osfc;
ssp126.end.flooded = ssp126.flooded(:,(end_start-2015)*12+1:(end_end-2015+1)*12);
ssp126.end.fh2osfc = ssp126.fh2osfc(:,(end_start-2015)*12+1:(end_end-2015+1)*12);
ssp126.end.swf     = ssp126.end.flooded + ssp126.end.fh2osfc;
ssp585.mid.flooded = ssp585.flooded(:,(mid_start-2015)*12+1:(mid_end-2015+1)*12);
ssp585.mid.fh2osfc = ssp585.fh2osfc(:,(mid_start-2015)*12+1:(mid_end-2015+1)*12);
ssp585.mid.swf     = ssp585.mid.flooded + ssp585.mid.fh2osfc;
ssp585.end.flooded = ssp585.flooded(:,(end_start-2015)*12+1:(end_end-2015+1)*12);
ssp585.end.fh2osfc = ssp585.fh2osfc(:,(end_start-2015)*12+1:(end_end-2015+1)*12);
ssp585.end.swf     = ssp585.end.flooded + ssp585.end.fh2osfc;

ssp126.bc.mid.swf(lakein,:)   = NaN; ssp126.bc.end.swf(lakein,:)   = NaN;
ssp585.bc.mid.swf(lakein,:)   = NaN; ssp585.bc.end.swf(lakein,:)   = NaN;
ssp126.mid.swf(lakein,:)      = NaN; ssp126.end.swf(lakein,:)      = NaN;
ssp585.mid.swf(lakein,:)      = NaN; ssp585.end.swf(lakein,:)      = NaN;
cmap = flipud(blue2red(17));

load('coastlines.mat');

[fig0,axs0,cb0,pos0] = plot_two_axes(xv,yv, ...
                    (nanmean(historical.bc.swf,2) - nanmean(historical.swf,2)).*area,...
                    (nanmean(historical.bc.swf,2) - nanmean(historical.swf,2)).*area,...
                    -21.25,21.25,'[km^2]');

[fig1,axs1,cb1,pos1] = plot_two_axes(xv,yv, ...
                    (nanmean(ssp126.bc.end.swf,2) - nanmean(ssp126.end.swf,2)).*area,...
                    (nanmean(ssp585.bc.end.swf,2) - nanmean(ssp585.end.swf,2)).*area,...
                    -21.25,21.25,'[km^2]');
t = add_title(axs1(1),'Impacts of BC on Surface Water',20,'out');
t.Position(2) = t.Position(2) + 0.025;

[fig2,axs2,cb2,pos2] = plot_two_axes(xv,yv, ...
                    (nanmean(ssp126.bc.end.swf,2) - nanmean(historical.bc.swf,2)).*area,...
                    (nanmean(ssp585.bc.end.swf,2) - nanmean(historical.bc.swf,2)).*area,...
                    -21.25,21.25,'[km^2]');

foc2d_ssp126 = NaN(720,360);
foc2d_ssp126(index_lnd) = (nanmean(ssp126.bc.end.swf,2) - nanmean(historical.bc.swf,2)).*area; % [km^2]
foc2d_ssp585 = NaN(720,360);
foc2d_ssp585(index_lnd) = (nanmean(ssp585.bc.end.swf,2) - nanmean(historical.bc.swf,2)).*area; % [km^2]
axs2(3) = axes('Position',[pos2(1) + pos2(3) + 0.02 pos2(2) 0.25 axs2(1).Position(4)+axs2(1).Position(2)-axs2(2).Position(2)]);

ctl        = NaN(720,360);
ssp126_end = NaN(720,360);
ssp585_end = NaN(720,360);
ctl(index_lnd) = nanmean(historical.bc.swf,2) .* area;
ssp126_end(index_lnd) = nanmean(ssp126.bc.end.swf,2) .* area;

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
                    (nanmean(ssp126.bc.end.swf,2) - nanmean(historical.swf,2)).*area,...
                    (nanmean(ssp585.bc.end.swf,2) - nanmean(historical.swf,2)).*area,...
                    -21.25,21.25,'[km^2]');

