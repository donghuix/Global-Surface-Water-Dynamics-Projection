clear;close all;clc;

re = 6.37122e6;% Earth radius

do_cal = false;

tag = 2; % tag = 1, only calibrate swf against glad
         % tag = 2, calibrate both swf and fsat. wt1 = 1, wt2 = 1
         % tag = 3, only calibrate fsat
         % tag = 4, calibrate both swf and fsat. wt1 = 1, wt2 = 3
         % tag = 5, calibrate both swf and fsat. wt1 = 1, wt2 = 2

if tag == 2
    wt1 = 1;
    wt2 = 1;
elseif tag == 4
    wt1 = 1;
    wt2 = 3;
elseif tag == 5
    wt1 = 1;
    wt2 = 2;
end

[e3sm_input, exportfig] = SetupEnvironment();

month_labels = {'01_Jan','02_Feb','03_Mar','04_Apr','05_May','06_Jun', ...
                '07_Jul','08_Aug','09_Sep','10_Oct','11_Nov','12_Dec'};

domain_file = [e3sm_input 'share/domains/domain.lnd.r05_oEC60to30v3.190418.nc'];
frac = ncread(domain_file,'frac'); frac = flipud(frac');
lon  = ncread(domain_file,'xc');
lat   = ncread(domain_file,'yc');
area2d= ncread(domain_file,'area');
cmap = getPanoply_cMap('EO_aura_omi_formal');
cmap = create_nonlinear_cmap(cmap,0,0.3,0.3,3);
load('coastlines.mat');
% Read GLAD seasonal surface water dynamics at 0.5x0.5

load('../data/GLAD05/GLAD05_permenant_water.mat');
load('index_lnd.mat');
xv = ncread('../data/domain_lnd_GLOBE_1d.nc','xv');
yv = ncread('../data/domain_lnd_GLOBE_1d.nc','yv');

% Read WAD2M, lon and lat is consistent with domain_file
Fw05 = ncread('../data/WAD2M_wetlands_2000-2020_05deg_Ver2.0.nc','Fw'); 
Fw05 = Fw05(:,:,(2001-2000)*12+1 : (2014-2000+1)*12);

pw(frac < 1) = NaN;
gladseason = NaN(length(index_lnd),12);
gladannual = NaN(720,360,12);
figure; set(gcf,'Position',[10 10 1200 800]);
for i = 1 : 12
    load(['../data/GLAD05/GLAD05_' month_labels{i} '.mat']);
    subplot_tight(3,4,i,[0.05 0.04]);

    GLAD05(frac < 1) = NaN;
    GLAD05 = GLAD05 - pw;
    GLAD05(GLAD05 <= 0) = 0;
    GLAD05 = fliplr(GLAD05');

    gladseason(:,i)   = GLAD05(index_lnd);
    gladannual(:,:,i) = GLAD05;
    patch(xv,yv,GLAD05(index_lnd)./100,'LineStyle','none'); clim([0 0.2]); hold on; grid on;
    ylim([-60 80]); set(gca,'XTick',[],'YTick',[]);
    plot(coastlon,coastlat,'k-','LineWidth',1);
    colormap(cmap);
end
save('../data/glad.mat','gladannual','gladseason');

gladannual = nanmean(gladannual,3);

Fwseason = NaN(length(index_lnd),12);
for i = 1 : 12
    tmp = nanmean(Fw05(:,:,i:12:end),3);
    Fwseason(:,i) = tmp(index_lnd);
end

% Load calibration simulations
numc = size(xv,2);
yr_start = 1973;
yr_end   = 2014;
yr_read  = 1993;
yr_glad  = 1999;

if exist(['../data/swf_cal_' num2str(tag) '.mat'],'file') && ~do_cal   
    load(['../data/swf_cal_' num2str(tag) '.mat']);
    load(['../data/par_cal_' num2str(tag) '.mat']);
    load(['../data/fsat_cal_' num2str(tag) '.mat'],'fsat_sea_cal');
else
    
    rng(4); % Random seed
    X = lhsdesign(100,2,'Criterion','correlation'); % fover, fc
    fover = exp(log(0.1) + (log(5)-log(0.1))*X(:,1));
    %fover = 0.1 + (0.5-0.1)*X(:,1);
    fc    = 0.001 + (0.4-0.001)*X(:,2);
%     figure;
%     plot(fc,fover,'bx','LineWidth',2); hold on;

    swf_sea_cal = NaN(numc,12);
    fsat_sea_cal= NaN(numc,12);

    swf_mon_cal = NaN(numc,(yr_end - yr_read+1)*12);
    fld_mon_cal = NaN(numc,(yr_end - yr_read+1)*12);
    h2o_mon_cal = NaN(numc,(yr_end - yr_read+1)*12);
    fsat_mon_cal= NaN(numc,(yr_end - yr_read+1)*12);
    
    swf_yr_cal = NaN(numc,yr_end - yr_read + 1);
    fld_yr_cal = NaN(numc,yr_end - yr_read + 1);
    h2o_yr_cal = NaN(numc,yr_end - yr_read + 1);
    fsat_yr_cal= NaN(numc,yr_end - yr_read + 1);
    
    NSEcal    = NaN(numc,1);
    fover_cal = ones(numc,1).*0.5;
    fc_cal    = ones(numc,1).*0.4;
    
    for i = 1 : 10
        load(['../data/cal_' num2str(i) '.mat']);
        swf = flooded + fh2osfc;
        if i < 10
            fprintf(['Calibration 0' num2str(i) ': ']);
        else
            fprintf(['Calibration '  num2str(i) ': ']);
        end
        for j = 1 : numc
            if mod(j,667*5) == 0
                fprintf([num2str(j/667) '%% ']);
            end
            if ~all(isnan(gladseason(j,:))) && ~all(gladseason(j,:) == 0)
                swfgrid = swf(j : numc : end, :);
                fldgrid = flooded(j : numc : end, :);
                h2ogrid = fh2osfc(j : numc : end, :);
                fsatgrid = fsat(j : numc :end, :);
    
                NSE = NaN(10,1);
                for k = 1 : 10
                    tmp = nanmean(reshape(swfgrid(k,(yr_glad - yr_read)*12+1:end)',[12,yr_end - yr_glad + 1]),2);
                    [~,~,NSE(k)] = estimate_evaluation_metric(gladseason(j,:)'./100,tmp(:));

                    tmp2 = nanmean(reshape(fsatgrid(k,(yr_glad - yr_read)*12+1:end)',[12,yr_end - yr_glad + 1]),2);
                    [~,~,NSE2(k)] = estimate_evaluation_metric(Fwseason(j,:)',tmp2(:));

                    if tag == 1
                        obj = NSE;
                    elseif tag == 2
                        obj = wt1.*NSE+wt2.*NSE2;
                    elseif tag == 3
                        obj = NSE2;
                    elseif tag == 4
                        obj = wt1.*NSE+wt2.*NSE2;
                    elseif tag == 5
                        obj = wt1.*NSE+wt2.*NSE2;
                    end
                end
                if i == 1
                    ind = find(obj == max(obj));
                    if ~isempty(ind)
                        ind = ind(1);

                        NSEcal(j) = obj(ind);
  
                        tmp  = reshape(swfgrid(ind,(yr_glad - yr_read)*12+1:end)',[12,yr_end - yr_glad + 1]);
                        swf_sea_cal(j,:) = nanmean(tmp,2);

                        tmp  = reshape(fsatgrid(ind,(yr_glad - yr_read)*12+1:end)',[12,yr_end - yr_glad + 1]);
                        fsat_sea_cal(j,:) = nanmean(tmp,2);

                        tmp  = reshape(swfgrid(ind,:)',[12,yr_end - yr_read + 1]);
                        swf_mon_cal(j,:) = swfgrid(ind,:);
                        swf_yr_cal(j,:)  = nanmean(tmp,1);

                        tmp1 = reshape(fldgrid(ind,:)',[12,yr_end - yr_read + 1]);
                        fld_mon_cal(j,:) = fldgrid(ind,:);
                        fld_yr_cal(j,:)  = nanmean(tmp1,1);

                        tmp2 = reshape(h2ogrid(ind,:)',[12,yr_end - yr_read + 1]);
                        h2o_mon_cal(j,:) = h2ogrid(ind,:);
                        h2o_yr_cal(j,:)  = nanmean(tmp2,1);

                        tmp3 = reshape(fsatgrid(ind,:)',[12,yr_end - yr_read + 1]);
                        fsat_mon_cal(j,:) = fsatgrid(ind,:);
                        fsat_yr_cal(j,:)  = nanmean(tmp3,1);

                        fover_cal(j) = fover(ind);
                        fc_cal(j)    = fc(ind);
                    end
                else
  
                    ind = find(obj == max(obj));

                    if ~isempty(ind)
                        if max(obj) > NSEcal(j)
                            ind = ind(1);
                            
                            NSEcal(j) = obj(ind);

                            tmp  = reshape(swfgrid(ind,(yr_glad - yr_read)*12+1:end)',[12,yr_end - yr_glad + 1]);
                            swf_sea_cal(j,:) = nanmean(tmp,2);
                            
                            tmp  = reshape(fsatgrid(ind,(yr_glad - yr_read)*12+1:end)',[12,yr_end - yr_glad + 1]);
                            fsat_sea_cal(j,:) = nanmean(tmp,2);

                            tmp = reshape(swfgrid(ind,:)',[12,yr_end - yr_read + 1]);
                            swf_mon_cal(j,:) = swfgrid(ind,:);
                            swf_yr_cal(j,:)  = nanmean(tmp,1);

                            tmp1 = reshape(fldgrid(ind,:)',[12,yr_end - yr_read + 1]);
                            fld_mon_cal(j,:) = fldgrid(ind,:);
                            fld_yr_cal(j,:)  = nanmean(tmp1,1);

                            tmp2 = reshape(h2ogrid(ind,:)',[12,yr_end - yr_read + 1]);
                            h2o_mon_cal(j,:) = h2ogrid(ind,:);
                            h2o_yr_cal(j,:)  = nanmean(tmp2,1);
                            
                            tmp3 = reshape(fsatgrid(ind,:)',[12,yr_end - yr_read + 1]);
                            fsat_mon_cal(j,:) = fsatgrid(ind,:);
                            fsat_yr_cal(j,:)  = nanmean(tmp3,1);

                            fover_cal(j) = fover((i-1)*10 + ind);
                            fc_cal(j)    = fc((i-1)*10 + ind);
                        end
                    end
                end
            end
        end
        fprintf('100%% \n');
    end
    
    save(['../data/swf_cal_' num2str(tag) '.mat'],'swf_mon_cal','swf_yr_cal', ...
                                                  'fld_mon_cal','fld_yr_cal', ...
                                                  'h2o_mon_cal','h2o_yr_cal', ...
                                                  'swf_sea_cal','-v7.3');
    save(['../data/fsat_cal_' num2str(tag) '.mat'],'fsat_mon_cal','fsat_yr_cal','fsat_sea_cal','-v7.3');
    save(['../data/par_cal_' num2str(tag) '.mat'],'fover_cal','fc_cal');

end

xc   = ncread('../data/domain_lnd_GLOBE_1d.nc','xc');
yc   = ncread('../data/domain_lnd_GLOBE_1d.nc','yc');
area = ncread('../data/domain_lnd_GLOBE_1d.nc','area');
load('LargeLakes.mat');
lakein = [];
for i = 1 : 20
    tmp = inpoly2([xc yc],[LargeLakes(i).X' LargeLakes(i).Y']);
    tmp = find(tmp == 1);
    lakein  = [lakein; tmp];
end

lakein2d = [];
for i = 1 : 20
    tmp = inpoly2([lon(:) lat(:)],[LargeLakes(i).X' LargeLakes(i).Y']);
    tmp = find(tmp == 1);
    lakein2d  = [lakein2d; tmp];
end

swf_sea_cal(lakein,:) = NaN;
fsat_sea_cal(lakein,:)= NaN;
swf_sea_cal(nanmean(swf_sea_cal,2) >= 0.3,:) = NaN;

S = shaperead('../data/TPBoundary_new(2021)/TPBoundary_new(2021).shp');

gladannual(lakein2d) = NaN;

figure; set(gcf,'Position',[10 10 1000 1200]);
axs(1) = subplot(2,1,1);
[cb,~,~] = plot_globalspatial(lon,lat,gladannual'./100,1,1); clim([0 0.3]); hold on; grid on;
colormap(gca,cmap);
ylim([-60 80]); xlim([-180 180]);
plot(S.X,S.Y,'k-','LineWidth',2);

pos1 = get(gca,'Position');

axs(2) = subplot(2,1,2);
swf_sea_cal(isnan(gladseason)) = NaN;
patch(xv,yv,nanmean(swf_sea_cal,2),'LineStyle','none'); clim([0 0.3]); hold on; grid on;
%plot(xc(nanmean(swfcal,2) > 0.3), yc(nanmean(swfcal,2) > 0.3),'rx','LineWidth',2);
ylim([-60 80]); xlim([-180 180]);%set(gca,'XTick',[],'YTick',[]); 
plot(coastlon,coastlat,'k-','LineWidth',1);
plot(S.X,S.Y,'k-','LineWidth',2);

colormap(gca,cmap);

axs(2).Position(2) = axs(2).Position(2) + 0.075;
axs(1).Position(1) = axs(1).Position(1) - 0.1;
axs(2).Position(1) = axs(2).Position(1) - 0.1;
axs(1).Position(3) = axs(1).Position(3) - 0.1;
axs(2).Position(3) = axs(2).Position(3) - 0.1;

pos2 = get(axs(2),'Position');
cb.Location = 'south';
cb.Position(1) = pos2(1);
cb.Position(2) = pos2(2)-0.05;
cb.Position(3) = pos2(3);
cb.Position(4) = 0.02;
cb.AxisLocation = 'out';
add_title(axs(1),'(a) GLAD',20,'out');
add_title(axs(2),'(b) Simulation',20,'out');

gladseason(lakein,:) = NaN;
[R2,RMSE,NSE,PBIAS] = estimate_evaluation_metric(nanmean(gladseason,2)./100.*area,nanmean(swf_sea_cal,2).*area);
KGE = estimateKGE(nanmean(gladseason,2)./100.*area,nanmean(swf_sea_cal,2).*area);

if tag ~= 3
str = {['\rho = ' num2str(round(sqrt(R2),2))] , ...
       ['NSE = ' num2str(round(NSE,2))]};
t = add_title(axs(2),str,18,'in');
t.Color = 'b';
t.Position(2) = t.Position(2) - 0.25;
end

axs(3) = axes('Position',[pos2(1) + pos2(3) + 0.02 pos2(2) 0.25 axs(1).Position(4)+axs(1).Position(2)-axs(2).Position(2)]);

sim = nanmean(swf_sea_cal,2).*area.*re^2./1e6;
obs = nanmean(gladseason,2)./100.*area.*re^2./1e6;
obs(isnan(sim)) = NaN;
sim2d = NaN(720,360);
obs2d = NaN(720,360);
sim2d(index_lnd) = sim;
obs2d(index_lnd) = obs;

plot(axs(3),nansum(obs2d,1),-89.75:0.5:89.75,'k-', 'LineWidth',2);hold on; grid on;
plot(axs(3),nansum(sim2d,1),-89.75:0.5:89.75,'r--','LineWidth',2); 

ylim([-60 80]);

set(axs(1),'Color',[0.8 0.8 0.8]);
set(axs(2),'Color',[0.8 0.8 0.8]);
leg = legend('GLAD','Simulation','FontSize',18,'FontWeight','bold');

[R2,RMSE,NSE,PBIAS,MSE,NSE1] = estimate_evaluation_metric(nansum(obs2d,1)',nansum(sim2d,1)');

str = {['\rho = ' num2str(round(sqrt(R2),2))] , ...
       ['NSE = ' num2str(round(NSE,2))]};
t = add_title(axs(3),str,18,'in');
t.Color = 'b';
t.Position(1) = t.Position(1) + 0.15;
t.Position(2) = t.Position(2) - 0.55;

figure; set(gcf,'Position',[10 10 1000 1200]);
axs(1) = subplot(2,1,1);
[cb,~,~] = plot_globalspatial(lon,lat,nanmean(Fw05,3)',1,1); clim([0 1]); hold on; grid on;
colormap(gca,cmap);
ylim([-60 80]); xlim([-180 180]);
plot(S.X,S.Y,'k-','LineWidth',2);

axs(2) = subplot(2,1,2);
% fsat_sea_cal(isnan(gladseason)) = NaN;
patch(xv,yv,nanmean(fsat_sea_cal,2),'LineStyle','none'); clim([0 1]); hold on; grid on;
colormap(gca,cmap);
ylim([-60 80]); xlim([-180 180]);%set(gca,'XTick',[],'YTick',[]); 
plot(coastlon,coastlat,'k-','LineWidth',1);
plot(S.X,S.Y,'k-','LineWidth',2);

axs(2).Position(2) = axs(2).Position(2) + 0.075;
axs(1).Position(1) = axs(1).Position(1) - 0.1;
axs(2).Position(1) = axs(2).Position(1) - 0.1;
axs(1).Position(3) = axs(1).Position(3) - 0.1;
axs(2).Position(3) = axs(2).Position(3) - 0.1;

pos2 = get(axs(2),'Position');
cb.Location = 'south';
cb.Position(1) = pos2(1);
cb.Position(2) = pos2(2)-0.05;
cb.Position(3) = pos2(3);
cb.Position(4) = 0.02;
cb.AxisLocation = 'out';
add_title(axs(1),'(a) WAD2M',20,'out');
add_title(axs(2),'(b) Simulation',20,'out');


axs(3) = axes('Position',[pos2(1) + pos2(3) + 0.02 pos2(2) 0.25 axs(1).Position(4)+axs(1).Position(2)-axs(2).Position(2)]);
fsat= nanmean(fsat_sea_cal,2).*area.*re^2./1e6;
fsat2d = NaN(720,360);
fw2d  = nanmean(Fw05,3).*area2d.*re^2./1e6;
fsat2d(index_lnd)= fsat;
fw2d(isnan(obs2d)) = NaN;

plot(axs(3),nansum(fw2d,1), -89.75:0.5:89.75,'k-','LineWidth',2); hold on;
plot(axs(3),nansum(fsat2d,1), -89.75:0.5:89.75,'b--','LineWidth',2); 
ylim([-60 80]);

set(axs(1),'Color',[0.8 0.8 0.8]);
set(axs(2),'Color',[0.8 0.8 0.8]);
leg = legend('WAD2M','Simulation','FontSize',18,'FontWeight','bold');

[R2,RMSE,NSE,PBIAS,MSE,NSE1] = estimate_evaluation_metric(nansum(fw2d,1)',nansum(fsat2d,1)');

str = {['\rho = ' num2str(round(sqrt(R2),2))] , ...
       ['NSE = ' num2str(round(NSE,2))]};
t = add_title(axs(3),str,18,'in');
t.Color = 'b';
t.Position(1) = t.Position(1) + 0.15;
t.Position(2) = t.Position(2) - 0.55;


% figure;
% fc_cal(isnan(nanmean(swf_sea_cal,2))) = NaN;
% fover_cal(isnan(nanmean(swf_sea_cal,2))) = NaN;
% subplot(2,1,1);
% patch(xv,yv,fover_cal,'LineStyle','none'); colorbar;
% subplot(2,1,2);
% patch(xv,yv,fc_cal,'LineStyle','none'); colorbar;
% 
% figure;
% plot(nanmean(gladseason,1)./100,'k-','LineWidth',2); hold on; grid on;
% plot(nanmean(swf_sea_cal,1),'r:','LineWidth',2);
% 
% figure;
% continent = struct([]);
% continent_code = {'af',    'ar',    'as',  'au',        'eu',    'gr',       'na',           'sa',           'si'     };
% continent_name = {'Africa','Arctic','Asia','Austrialia','Europe','Greenland','North America','South America','Siberia'};
% 
% subplot(3,3,1);
% plot(nanmean(gladseason,1)./100,'k-','LineWidth',2); hold on; grid on;
% plot(nanmean(swf_sea_cal,1),'r--','LineWidth',2);
% add_title(gca,'Global');
% 
% k = 2;
% for i = [1 2 3 4 5 7 8 9]
%    code = continent_code{i};
%    continent(i).code  = code;
%    continent(i).name  = continent_name{i};
%    continent(i).index = []; 
%    continent(i).swf   = NaN(20,1);
%    S = shaperead(['../data/HydroBASINS/hybas_' code '_lev01-06_v1c/hybas_' code '_lev01_v1c.shp']);
%    for j = 1 : length(S)
%        tmp = inpoly2([xc(:) yc(:)],[S(j).X' S(j).Y']);
%        tmp = find(tmp == 1);
%        %plot(S(j).X,S(j).Y,'-','LineWidth',2);
%        continent(i).index = [continent(i).index; tmp];
%    end
% 
%    subplot(3,3,k);
%    plot(nanmean(gladseason(continent(i).index,:),1)./100,'k-','LineWidth',2); hold on; grid on;
%    plot(nanmean(swf_sea_cal(continent(i).index,:),1),'r--','LineWidth',2);
% 
%    add_title(gca,continent_name{i});
% 
%    k = k + 1;
% end
% 
% figure;
% for i = 1 : length(continent_code)
%     subplot(3,3,i);
%     plot(nanmean(swf_yr_cal(continent(i).index,:),1),'r--','LineWidth',2);
%     add_title(gca,continent_name{i});
% end
% 
% figure;
% patch(xv,yv,nanmean(swf_sea_cal,2) - nanmean(gladseason./100,2),'LineStyle','none'); 
% clim([-0.1 0.1]); hold on; grid on;
% colormap(blue2red(11));

% for i = 1 : 50
%     if i <= 10
%         plot(LargeLakes(i).X,LargeLakes(i).Y,'r-','LineWidth',2);
%     elseif i <= 20
%         plot(LargeLakes(i).X,LargeLakes(i).Y,'b-','LineWidth',2);
%     elseif i <= 30
%         plot(LargeLakes(i).X,LargeLakes(i).Y,'g-','LineWidth',2);
%     elseif i <= 40
%         plot(LargeLakes(i).X,LargeLakes(i).Y,'m-','LineWidth',2);
%     else
%         plot(LargeLakes(i).X,LargeLakes(i).Y,'c-','LineWidth',2);
%     end
% end
