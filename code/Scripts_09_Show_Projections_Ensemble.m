clear;close all;clc;

[e3sm_input, exportfig] = SetupEnvironment();

models = {'gfdl-esm4', 'ipsl-cm6a-lr', 'mpi-esm1-2-hr', 'mri-esm2-0', 'ukesm1-0-ll'};
scenarios = {'ssp126','ssp585'};
times  = {'mid','end'};
mid_start = 2041;
mid_end   = 2070;
end_start = 2071;
end_end   = 2100;
month_labels = {'01_Jan','02_Feb','03_Mar','04_Apr','05_May','06_Jun', ...
                '07_Jul','08_Aug','09_Sep','10_Oct','11_Nov','12_Dec'};

load('coastlines.mat');

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

if exist('ensemble.mat','file')
    load('ensemble.mat');
else
    wetland_ensemble = struct([]);
    global_ensemble  = struct([]);
    continent = struct([]);
    continent_code = {'af',    'ar',    'as',  'au',        'eu',    'gr',       'na',           'sa',           'si'     };
    continent_name = {'Africa','Arctic','Asia','Austrialia','Europe','Greenland','North America','South America','Siberia'};
    
    k = 1;
    for i = [1 2 3 4 5 7 8 9]
       code = continent_code{i};
       continent(k).code  = code;
       continent(k).name  = continent_name{i};
       continent(k).index = []; 
       S = shaperead(['../data/HydroBASINS/hybas_' code '_lev01-06_v1c/hybas_' code '_lev01_v1c.shp']);
       for j = 1 : length(S)
           tmp = inpoly2([xc(:) yc(:)],[S(j).X' S(j).Y']);
           tmp = find(tmp == 1);
           %plot(S(j).X,S(j).Y,'-','LineWidth',2);
           continent(k).index = [continent(k).index; tmp];
       end
       k = k + 1;
    end

    load(['../data/fsat_cal_12.mat'],'fsat_sea_cal');
    irm = find(isnan(nanmean(fsat_sea_cal,2)));
    for i = 1 : length(models)
        model = models{i};
        
        disp(model);
        historical = load(['../projection/projection_cal12_' model '_historical.mat'],'flooded','fh2osfc','fsat','tsa','rain');
        ssp126     = load(['../projection/projection_cal12_' model '_ssp126.mat'],'flooded','fh2osfc','fsat','tsa','rain');
        ssp585     = load(['../projection/projection_cal12_' model '_ssp585.mat'],'flooded','fh2osfc','fsat','tsa','rain');
        global_ensemble(1).CTL.swf(:,i)   = nansum((historical.flooded+historical.fh2osfc).*area);
        global_ensemble(1).CTL.fsat(:,i)  = nansum(historical.fsat.*area);
        for ii = 1 : 8
            continent(ii).CTL.swf(:,i)    = nansum((historical.flooded(continent(ii).index,:) +  ...
                                                    historical.fh2osfc(continent(ii).index,:)).* ...
                                                    area(continent(ii).index));
            continent(ii).CTL.fsat(:,i)   = nansum(historical.fsat(continent(ii).index,:).*area(continent(ii).index));
            continent(ii).CTL.tsa(:,i)    = nansum(historical.tsa(continent(ii).index,:).*area(continent(ii).index))./ ...
                                                   nansum(area(continent(ii).index));
        end

        historical.flooded = historical.flooded(:,1:360);
        historical.fh2osfc = historical.fh2osfc(:,1:360);
        historical.swf     = historical.flooded + historical.fh2osfc;
        historical.fsat    = historical.fsat(:,1:360);
        historical.tsa     = historical.tsa(:,1:360);
        historical.rain     = historical.rain(:,1:360);
        wetland_ensemble(1).CTL.swf(:,i)  = nanmean(historical.swf,2);
        wetland_ensemble(1).CTL.fsat(:,i) = nanmean(historical.fsat,2);
        wetland_ensemble(1).CTL.tsa(:,i)  = nanmean(historical.tsa,2);
        wetland_ensemble(1).CTL.rain(:,i) = nanmean(historical.rain,2);
        
    
        flooded = ssp126.flooded(:,(mid_start-2015)*12+1:(mid_end-2015+1)*12);
        fh2osfc = ssp126.fh2osfc(:,(mid_start-2015)*12+1:(mid_end-2015+1)*12);
        fsat    = ssp126.fsat(:,(mid_start-2015)*12+1:(mid_end-2015+1)*12);
        tsa     = ssp126.tsa(:,(mid_start-2015)*12+1:(mid_end-2015+1)*12);
        rain    = ssp126.rain(:,(mid_start-2015)*12+1:(mid_end-2015+1)*12);
        wetland_ensemble(1).ssp126.mid.swf(:,i)  = nanmean(flooded + fh2osfc,2);
        wetland_ensemble(1).ssp126.mid.fsat(:,i) = nanmean(fsat,2);
        wetland_ensemble(1).ssp126.mid.tsa(:,i)  = nanmean(tsa,2);
        wetland_ensemble(1).ssp126.mid.rain(:,i) = nanmean(rain,2);
        global_ensemble(1).ssp126.swf(:,i)   = nansum((ssp126.flooded+ssp126.fh2osfc).*area);
        global_ensemble(1).ssp126.fsat(:,i)  = nansum(ssp126.fsat.*area);
        for ii = 1 : 8
            continent(ii).ssp126.swf(:,i)    = nansum((ssp126.flooded(continent(ii).index,:) +  ...
                                                    ssp126.fh2osfc(continent(ii).index,:)).* ...
                                                    area(continent(ii).index));
            continent(ii).ssp126.fsat(:,i)   = nansum(ssp126.fsat(continent(ii).index,:).*area(continent(ii).index));
            continent(ii).ssp126.tsa(:,i)    = nansum(ssp126.tsa(continent(ii).index,:).*area(continent(ii).index))./ ...
                                                       nansum(area(continent(ii).index));
        end

        flooded = ssp126.flooded(:,(end_start-2015)*12+1:(end_end-2015+1)*12);
        fh2osfc = ssp126.fh2osfc(:,(end_start-2015)*12+1:(end_end-2015+1)*12);
        fsat    = ssp126.fsat(:,(end_start-2015)*12+1:(end_end-2015+1)*12);
        tsa     = ssp126.tsa(:,(end_start-2015)*12+1:(end_end-2015+1)*12);
        rain    = ssp126.rain(:,(end_start-2015)*12+1:(end_end-2015+1)*12);
        wetland_ensemble(1).ssp126.end.swf(:,i)  = nanmean(flooded + fh2osfc,2);
        wetland_ensemble(1).ssp126.end.fsat(:,i) = nanmean(fsat,2);
        wetland_ensemble(1).ssp126.end.tsa(:,i)  = nanmean(tsa,2);
        wetland_ensemble(1).ssp126.end.rain(:,i) = nanmean(rain,2);
    
        flooded = ssp585.flooded(:,(mid_start-2015)*12+1:(mid_end-2015+1)*12);
        fh2osfc = ssp585.fh2osfc(:,(mid_start-2015)*12+1:(mid_end-2015+1)*12);
        fsat    = ssp585.fsat(:,(mid_start-2015)*12+1:(mid_end-2015+1)*12);
        tsa     = ssp585.tsa(:,(mid_start-2015)*12+1:(mid_end-2015+1)*12);
        rain    = ssp585.rain(:,(mid_start-2015)*12+1:(mid_end-2015+1)*12);
        wetland_ensemble(1).ssp585.mid.swf(:,i)  = nanmean(flooded + fh2osfc,2);
        wetland_ensemble(1).ssp585.mid.fsat(:,i) = nanmean(fsat,2);
        wetland_ensemble(1).ssp585.mid.tsa(:,i)  = nanmean(tsa,2);
        wetland_ensemble(1).ssp585.mid.rain(:,i) = nanmean(rain,2);
        global_ensemble(1).ssp585.swf(:,i)   = nansum((ssp585.flooded+ssp585.fh2osfc).*area);
        global_ensemble(1).ssp585.fsat(:,i)  = nansum(ssp585.fsat.*area);
        for ii = 1 : 8
            continent(ii).ssp585.swf(:,i)    = nansum((ssp585.flooded(continent(ii).index,:) +  ...
                                                    ssp585.fh2osfc(continent(ii).index,:)).* ...
                                                    area(continent(ii).index));
            continent(ii).ssp585.fsat(:,i)   = nansum(ssp585.fsat(continent(ii).index,:).*area(continent(ii).index));
            continent(ii).ssp585.tsa(:,i)    = nansum(ssp585.tsa(continent(ii).index,:).*area(continent(ii).index))./ ...
                                                       nansum(area(continent(ii).index));
        end

        flooded = ssp585.flooded(:,(end_start-2015)*12+1:(end_end-2015+1)*12);
        fh2osfc = ssp585.fh2osfc(:,(end_start-2015)*12+1:(end_end-2015+1)*12);
        fsat    = ssp585.fsat(:,(end_start-2015)*12+1:(end_end-2015+1)*12);
        tsa     = ssp585.tsa(:,(end_start-2015)*12+1:(end_end-2015+1)*12);
        rain    = ssp585.rain(:,(end_start-2015)*12+1:(end_end-2015+1)*12);
        wetland_ensemble(1).ssp585.end.swf(:,i)  = nanmean(flooded + fh2osfc,2);
        wetland_ensemble(1).ssp585.end.fsat(:,i) = nanmean(fsat,2);
        wetland_ensemble(1).ssp585.end.tsa(:,i)  = nanmean(tsa,2);
        wetland_ensemble(1).ssp585.end.rain(:,i) = nanmean(rain,2);
    
        wetland_ensemble(1).ssp126.mid.swf(lakein,i)   = NaN; 
        wetland_ensemble(1).ssp126.end.swf(lakein,i)   = NaN;
        wetland_ensemble(1).ssp126.mid.fsat(lakein,i)  = NaN; 
        wetland_ensemble(1).ssp126.end.fsat(lakein,i)  = NaN;
        wetland_ensemble(1).ssp585.mid.swf(lakein,i)   = NaN; 
        wetland_ensemble(1).ssp585.end.swf(lakein,i)   = NaN;
        wetland_ensemble(1).ssp585.mid.fsat(lakein,i)  = NaN; 
        wetland_ensemble(1).ssp585.end.fsat(lakein,i)  = NaN;
        
        wetland_ensemble(1).CTL.fsat(irm,i) = 0;
        wetland_ensemble(1).ssp126.mid.fsat(irm,i) = 0;
        wetland_ensemble(1).ssp126.end.fsat(irm,i) = 0;
        wetland_ensemble(1).ssp585.mid.fsat(irm,i) = 0;
        wetland_ensemble(1).ssp585.end.fsat(irm,i) = 0;
    
    end
    
    save('ensemble.mat','global_ensemble','wetland_ensemble','continent');
end

% FSAT
[fig1,axs1,cb1,pos1] = plot_two_axes(xv,yv, ...
                    nanmean(wetland_ensemble.ssp126.end.fsat-wetland_ensemble.CTL.fsat,2).*area,...
                    nanmean(wetland_ensemble.ssp585.end.fsat-wetland_ensemble.CTL.fsat,2).*area,...
                    -100, 100,'[km^2]');
colormap(axs1(1),cmap_new);
colormap(axs1(2),cmap_new);
add_title(axs1(1),'SSP126',20,'out');
add_title(axs1(2),'SSP585',20,'out');

axs1(3) = axes('Position',[pos1(1) + pos1(3) + 0.02 pos1(2) 0.25 axs1(1).Position(4)+axs1(1).Position(2)-axs1(2).Position(2)]);

ctl2d_ens        = NaN(360,5);
foc2d_ssp126_ens = NaN(360,5);
foc2d_ssp585_ens = NaN(360,5);
for i = 1 : 5
    
    ctl2d        = NaN(720,360);
    foc2d_ssp126 = NaN(720,360);
    foc2d_ssp585 = NaN(720,360);
    ctl2d(index_lnd)        =  wetland_ensemble.CTL.fsat(:,i).*area;
    foc2d_ssp126(index_lnd) = (wetland_ensemble.ssp126.end.fsat(:,i) - wetland_ensemble.CTL.fsat(:,i)).*area; % [km^2]
    foc2d_ssp585(index_lnd) = (wetland_ensemble.ssp585.end.fsat(:,i) - wetland_ensemble.CTL.fsat(:,i)).*area; % [km^2]
    
    ctl2d(isnan(gladannual))        = NaN;
    foc2d_ssp126(isnan(gladannual)) = NaN;
    foc2d_ssp585(isnan(gladannual)) = NaN;
    
    ctl2d_ens(:,i)        = movmean(nansum(ctl2d       ,1),5);
    foc2d_ssp126_ens(:,i) = movmean(nansum(foc2d_ssp126,1),5);
    foc2d_ssp585_ens(:,i) = movmean(nansum(foc2d_ssp585,1),5);
    %plot(axs1(3),movmean(nansum(foc2d_ssp126,1),5)./movmean(nansum(ctl2d,1),5),-89.75:0.5:89.75,'-', 'LineWidth',2);hold on; grid on;    %plot(axs1(3),zeros(360,1),-89.75:0.5:89.75,'k--','LineWidth',1); 
end

x = -89.75:0.5:89.75;
y1 = max(foc2d_ssp126_ens./ctl2d_ens,[],2)';
y2 = min(foc2d_ssp126_ens./ctl2d_ens,[],2)';
p = patch(axs1(3),[y1(69:332) fliplr(y2(69:332))],[x(69:332) fliplr(x(69:332))],'r'); hold on; grid on;
p.FaceColor = [49 67 143]./255; p.FaceAlpha = 0.3; p.EdgeColor = 'none';

y1 = max(foc2d_ssp585_ens./ctl2d_ens,[],2)';
y2 = min(foc2d_ssp585_ens./ctl2d_ens,[],2)';
p = patch(axs1(3),[y1(69:332) fliplr(y2(69:332))],[x(69:332) fliplr(x(69:332))],'b'); hold on; grid on;
p.FaceColor = [194 62 94]./255; p.FaceAlpha = 0.3; p.EdgeColor = 'none';
h(1) = plot(axs1(3),mean(foc2d_ssp126_ens(69:332,:)./ctl2d_ens(69:332,:),2),x(69:332),'-','Color',[49 67 143]./255,'LineWidth',3);
h(2) = plot(axs1(3),mean(foc2d_ssp585_ens(69:332,:)./ctl2d_ens(69:332,:),2),x(69:332),'-','Color',[194 62 94]./255,'LineWidth',3);
legend(h,{'SSP126','SSP585'},'FontSize',15,'FontWeight','bold');
set(axs1(3),'FontSize',13);xlim(axs1(3),[-0.6 0.6]); xticks(axs1(3),[-0.6 : 0.2 : 0.6]);
% %plot(axs1(3),zeros(360,1),-89.75:0.5:89.75,'k--','LineWidth',1); 
ylim(axs1(3),[x(69) x(332)]);
ylim(axs1(1),[x(69) x(332)]);
ylim(axs1(2),[x(69) x(332)]);


% SWF
[fig2,axs2,cb2,pos2] = plot_two_axes(xv,yv, ...
                    nanmean(wetland_ensemble.ssp126.end.swf-wetland_ensemble.CTL.swf,2).*area,...
                    nanmean(wetland_ensemble.ssp585.end.swf-wetland_ensemble.CTL.swf,2).*area,...
                    -30, 30,'[km^2]');
colormap(axs2(1),cmap_new);
colormap(axs2(2),cmap_new);
add_title(axs2(1),'SSP126',20,'out');
add_title(axs2(2),'SSP585',20,'out');

axs2(3) = axes('Position',[pos2(1) + pos2(3) + 0.02 pos2(2) 0.25 axs2(1).Position(4)+axs2(1).Position(2)-axs2(2).Position(2)]);

ctl2d_ens        = NaN(360,5);
foc2d_ssp126_ens = NaN(360,5);
foc2d_ssp585_ens = NaN(360,5);
for i = 1 : 5
    
    ctl2d        = NaN(720,360);
    foc2d_ssp126 = NaN(720,360);
    foc2d_ssp585 = NaN(720,360);
    ctl2d(index_lnd)        =  wetland_ensemble.CTL.swf(:,i).*area;
    foc2d_ssp126(index_lnd) = (wetland_ensemble.ssp126.end.swf(:,i) - wetland_ensemble.CTL.swf(:,i)).*area; % [km^2]
    foc2d_ssp585(index_lnd) = (wetland_ensemble.ssp585.end.swf(:,i) - wetland_ensemble.CTL.swf(:,i)).*area; % [km^2]
    
    ctl2d(isnan(gladannual))        = NaN;
    foc2d_ssp126(isnan(gladannual)) = NaN;
    foc2d_ssp585(isnan(gladannual)) = NaN;
    
    ctl2d_ens(:,i)        = movmean(nansum(ctl2d       ,1),5);
    foc2d_ssp126_ens(:,i) = movmean(nansum(foc2d_ssp126,1),5);
    foc2d_ssp585_ens(:,i) = movmean(nansum(foc2d_ssp585,1),5);
    %plot(axs1(3),movmean(nansum(foc2d_ssp126,1),5)./movmean(nansum(ctl2d,1),5),-89.75:0.5:89.75,'-', 'LineWidth',2);hold on; grid on;    %plot(axs1(3),zeros(360,1),-89.75:0.5:89.75,'k--','LineWidth',1); 
end

x = -89.75:0.5:89.75;
y1 = max(foc2d_ssp126_ens./ctl2d_ens,[],2)';
y2 = min(foc2d_ssp126_ens./ctl2d_ens,[],2)';
p = patch(axs2(3),[y1(69:332) fliplr(y2(69:332))],[x(69:332) fliplr(x(69:332))],'r'); hold on; grid on;
p.FaceColor = [49 67 143]./255; p.FaceAlpha = 0.3; p.EdgeColor = 'none';

y1 = max(foc2d_ssp585_ens./ctl2d_ens,[],2)';
y2 = min(foc2d_ssp585_ens./ctl2d_ens,[],2)';
p = patch(axs2(3),[y1(69:332) fliplr(y2(69:332))],[x(69:332) fliplr(x(69:332))],'b'); hold on; grid on;
p.FaceColor = [194 62 94]./255; p.FaceAlpha = 0.3; p.EdgeColor = 'none';
h(1) = plot(axs2(3),mean(foc2d_ssp126_ens(69:332,:),2),x(69:332),'-','Color',[49 67 143]./255,'LineWidth',3);
h(2) = plot(axs2(3),mean(foc2d_ssp585_ens(69:332,:),2),x(69:332),'-','Color',[194 62 94]./255,'LineWidth',3);
legend(h,{'SSP126','SSP585'},'FontSize',15,'FontWeight','bold');
set(axs2(3),'FontSize',13);%xlim(axs2(3),[-1 1]); xticks(axs2(3),[-1 : 0.2 : 1]);
% %plot(axs1(3),zeros(360,1),-89.75:0.5:89.75,'k--','LineWidth',1); 
ylim(axs2(3),[x(69) x(332)]);
ylim(axs2(1),[x(69) x(332)]);
ylim(axs2(2),[x(69) x(332)]);

figure;
axs(1) = subplot(2,1,1);

a = nanmean(reshape(nanmean(global_ensemble.CTL.fsat,2),[12 44]),1)./1e6;
b = nanmean(reshape(min(global_ensemble.CTL.fsat,[],2),[12 44]),1)./1e6;
c = nanmean(reshape(max(global_ensemble.CTL.fsat,[],2),[12 44]),1)./1e6;
x = 1971 : 2014;
patch([x'; flip(x')], [b'; flip(c')], 'k', 'FaceAlpha',0.2,'EdgeColor','none'); hold on; grid on;
plot(x,a,'k-','LineWidth',2);

a = nanmean(reshape(nanmean(global_ensemble.ssp126.fsat,2),[12 86]),1)./1e6;
b = nanmean(reshape(min(global_ensemble.ssp126.fsat,[],2),[12 86]),1)./1e6;
c = nanmean(reshape(max(global_ensemble.ssp126.fsat,[],2),[12 86]),1)./1e6;
x = 2015 : 2100;
patch([x'; flip(x')], [b'; flip(c')], 'b', 'FaceAlpha',0.2,'EdgeColor','none'); hold on; grid on;
plot(x,a,'b-','LineWidth',2);

a = nanmean(reshape(nanmean(global_ensemble.ssp585.fsat,2),[12 86]),1)./1e6;
b = nanmean(reshape(min(global_ensemble.ssp585.fsat,[],2),[12 86]),1)./1e6;
c = nanmean(reshape(max(global_ensemble.ssp585.fsat,[],2),[12 86]),1)./1e6;
x = 2015 : 2100;
patch([x'; flip(x')], [b'; flip(c')], 'r', 'FaceAlpha',0.2,'EdgeColor','none'); hold on;
plot(x,a,'r-','LineWidth',2);
xlim([1971 2100]);
xticklabels({''})

axs(2) = subplot(2,1,2);

a = nanmean(reshape(nanmean(global_ensemble.CTL.swf,2),[12 44]),1)./1e6;
b = nanmean(reshape(min(global_ensemble.CTL.swf,[],2),[12 44]),1)./1e6;
c = nanmean(reshape(max(global_ensemble.CTL.swf,[],2),[12 44]),1)./1e6;
x = 1971 : 2014;
patch([x'; flip(x')], [b'; flip(c')], 'k', 'FaceAlpha',0.2,'EdgeColor','none'); hold on; grid on;
plot(x,a,'k-','LineWidth',2);

a = nanmean(reshape(nanmean(global_ensemble.ssp126.swf,2),[12 86]),1)./1e6;
b = nanmean(reshape(min(global_ensemble.ssp126.swf,[],2),[12 86]),1)./1e6;
c = nanmean(reshape(max(global_ensemble.ssp126.swf,[],2),[12 86]),1)./1e6;
x = 2015 : 2100;
patch([x'; flip(x')], [b'; flip(c')], 'b', 'FaceAlpha',0.2,'EdgeColor','none'); hold on; grid on;
plot(x,a,'b-','LineWidth',2);

a = nanmean(reshape(nanmean(global_ensemble.ssp585.swf,2),[12 86]),1)./1e6;
b = nanmean(reshape(min(global_ensemble.ssp585.swf,[],2),[12 86]),1)./1e6;
c = nanmean(reshape(max(global_ensemble.ssp585.swf,[],2),[12 86]),1)./1e6;
x = 2015 : 2100;
patch([x'; flip(x')], [b'; flip(c')], 'r', 'FaceAlpha',0.2,'EdgeColor','none'); hold on;
plot(x,a,'r-','LineWidth',2);
%xlim([x(1) x(end)]);
xlim([1971 2100]);
axs(1).Position(2) = axs(1).Position(2) - 0.05;
add_title(axs(1),'(a). Total Wetland',20,'out');
%add_title(axs(2),'(c). Saturated Wetland',20,'out');
add_title(axs(2),'(b). Inundated Wetland',20,'out');

set(axs(1),'FontSize',15); set(axs(2),'FontSize',15); %set(axs(3),'FontSize',15); 

han=axes(gcf,'visible','off'); 
han.YLabel.Visible='on';
ylabel(han,'million km^{2}','FontSize',18,'FontWeight','bold');
han.Position(1) = han.Position(1) - 0.025;

figure;
for i = 1 : 8
    subplot(2,4,i);
    a = nanmean(reshape(nanmean(continent(i).CTL.swf,2),[12 44]),1)./1e6;
    b = nanmean(reshape(min(continent(i).CTL.swf,[],2),[12 44]),1)./1e6;
    c = nanmean(reshape(max(continent(i).CTL.swf,[],2),[12 44]),1)./1e6;
    x = 1971 : 2014;
    patch([x'; flip(x')], [b'; flip(c')], 'k', 'FaceAlpha',0.2,'EdgeColor','none'); hold on; grid on;
    plot(x,a,'k-','LineWidth',2);
    
    a = nanmean(reshape(nanmean(continent(i).ssp126.swf,2),[12 86]),1)./1e6;
    b = nanmean(reshape(min(continent(i).ssp126.swf,[],2),[12 86]),1)./1e6;
    c = nanmean(reshape(max(continent(i).ssp126.swf,[],2),[12 86]),1)./1e6;
    x = 2015 : 2100;
    patch([x'; flip(x')], [b'; flip(c')], 'b', 'FaceAlpha',0.2,'EdgeColor','none'); hold on; grid on;
    plot(x,a,'b-','LineWidth',2);
    
    a = nanmean(reshape(nanmean(continent(i).ssp585.swf,2),[12 86]),1)./1e6;
    b = nanmean(reshape(min(continent(i).ssp585.swf,[],2),[12 86]),1)./1e6;
    c = nanmean(reshape(max(continent(i).ssp585.swf,[],2),[12 86]),1)./1e6;
    x = 2015 : 2100;
    patch([x'; flip(x')], [b'; flip(c')], 'r', 'FaceAlpha',0.2,'EdgeColor','none'); hold on;
    plot(x,a,'r-','LineWidth',2);
    %xlim([x(1) x(end)]);
    xlim([1971 2100]);
    add_title(gca,continent(i).name,15,'out');
end

figure;
for i = 1 : 8
    subplot(2,4,i);
    a = nanmean(reshape(nanmean(continent(i).CTL.tsa,2),[12 44]),1)-273.15;
    b = nanmean(reshape(min(continent(i).CTL.tsa,[],2),[12 44]),1)-273.15;
    c = nanmean(reshape(max(continent(i).CTL.tsa,[],2),[12 44]),1)-273.15;
    x = 1971 : 2014;
    patch([x'; flip(x')], [b'; flip(c')], 'k', 'FaceAlpha',0.2,'EdgeColor','none'); hold on; grid on;
    plot(x,a,'k-','LineWidth',2);
    
    a = nanmean(reshape(nanmean(continent(i).ssp126.tsa,2),[12 86]),1)-273.15;
    b = nanmean(reshape(min(continent(i).ssp126.tsa,[],2),[12 86]),1)-273.15;
    c = nanmean(reshape(max(continent(i).ssp126.tsa,[],2),[12 86]),1)-273.15;
    x = 2015 : 2100;
    patch([x'; flip(x')], [b'; flip(c')], 'b', 'FaceAlpha',0.2,'EdgeColor','none'); hold on; grid on;
    plot(x,a,'b-','LineWidth',2);
    
    a = nanmean(reshape(nanmean(continent(i).ssp585.tsa,2),[12 86]),1)-273.15;
    b = nanmean(reshape(min(continent(i).ssp585.tsa,[],2),[12 86]),1)-273.15;
    c = nanmean(reshape(max(continent(i).ssp585.tsa,[],2),[12 86]),1)-273.15;
    x = 2015 : 2100;
    patch([x'; flip(x')], [b'; flip(c')], 'r', 'FaceAlpha',0.2,'EdgeColor','none'); hold on;
    plot(x,a,'r-','LineWidth',2);
    %xlim([x(1) x(end)]);
    xlim([1971 2100]);
    add_title(gca,continent(i).name,15,'out');
end
% pc = ncread('/Users/xudo627/Downloads/gpw-v4-population-count-rev11_totpop_30_min_nc/gpw_v4_population_count_rev11_30_min.nc','Population Count, v4.11 (2000, 2005, 2010, 2015, 2020): 30 arc-minutes');
% pc = fliplr(pc(:,:,5));
% 
% pc  = pc(index_lnd);
% foc1 = nanmean(wetland_ensemble.ssp126.end.swf-wetland_ensemble.CTL.swf,2).*area;
% foc2 = nanmean(wetland_ensemble.ssp585.end.swf-wetland_ensemble.CTL.swf,2).*area;
% 
% figure;
% pc(pc <= 10) = NaN;
% patch(xv,yv,nanmean(wetland_ensemble.CTL.swf,2).*area./pc.*1e6,'LineStyle','none'); hold on; colorbar; 
% set(gca,'ColorScale','log');clim([1e0 1e7]);
% 
% 
% figure;
% pc(pc <= 10) = NaN;
% patch(xv,yv,nanmean(wetland_ensemble.ssp585.end.swf-wetland_ensemble.CTL.swf,2).*area./pc.*1e6,'LineStyle','none'); hold on; colorbar; 
% colormap(flipud(blue2red(11))); clim([-1e7 1e7]);
% 
% a = nanmean(wetland_ensemble.CTL.swf,2).*area./pc.*1e6;
% b = nanmean(wetland_ensemble.ssp585.end.swf,2).*area./pc.*1e6;
% 
% figure;
% ecdf(a); hold on;
% ecdf(b);


