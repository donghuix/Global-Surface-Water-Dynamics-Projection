clear;close all;clc;

% Script is to plot the validation of calibrated model

re  = 6.37122e6;% Earth radius
tag = 12;
[e3sm_input, exportfig] = SetupEnvironment();

month_labels = {'01_Jan','02_Feb','03_Mar','04_Apr','05_May','06_Jun', ...
                '07_Jul','08_Aug','09_Sep','10_Oct','11_Nov','12_Dec'};

labels = {'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)','(j)'};

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
end
gladannual = nanmean(gladannual,3);

load(['../data/swf_cal_' num2str(tag) '.mat']);

xc   = ncread('../data/domain_lnd_GLOBE_1d.nc','xc');
yc   = ncread('../data/domain_lnd_GLOBE_1d.nc','yc');
area = ncread('../data/domain_lnd_GLOBE_1d.nc','area');
area = area.*(re^2) ./ 1e6; % square km
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
gladseason(lakein,:) = NaN;
fsat_sea_cal(lakein,:)= NaN;
gladannual(lakein2d) = NaN;

ind = find(nanmean(swf_sea_cal,2) > 0.3 & nanmean(gladseason./100,2) < 0.3);
save('problem.mat','ind');

d = nanmean(swf_sea_cal,2) - nanmean(gladseason./100,2);
ind2 = find(d >= 0.25);
%swf_sea_cal(d >= 0.25,:) = NaN;

figure; set(gcf,'Position',[10 10 1000 1200]);
axs(1) = subplot(2,1,1);
swf_sea_cal(isnan(gladseason)) = NaN;
patch(xv,yv,nanmean(swf_sea_cal,2),'LineStyle','none'); clim([0 0.3]); hold on; grid on;
%plot(xc(ind),yc(ind),'rx','LineWidth',2);
%plot(xc(nanmean(gladseason./100,2) > 0.3),yc(nanmean(gladseason./100,2) > 0.3),'gx','LineWidth',2);
%plot(xc(nanmean(swfcal,2) > 0.3), yc(nanmean(swfcal,2) > 0.3),'rx','LineWidth',2);
ylim([-60 80]); xlim([-180 180]);%set(gca,'XTick',[],'YTick',[]); 
plot(coastlon,coastlat,'k-','LineWidth',1);

colormap(gca,cmap); cb = colorbar;

axs(1).Position(1) = axs(1).Position(1) - 0.05;
axs(1).Position(3) = axs(1).Position(3) - 0.1;

pos1 = get(axs(1),'Position');
cb.Location = 'south';
cb.Position(1) = pos1(1);
cb.Position(2) = pos1(2)-0.05;
cb.Position(3) = pos1(3);
cb.Position(4) = 0.02;
cb.AxisLocation = 'out'; cb.FontSize = 15;

add_title(axs(1),'(a) Simulation',20,'out');

gladseason(lakein,:) = NaN;
swf_sea_cal(nanmean(swf_sea_cal,2) >= 0.3 & nanmean(gladseason./100,2) < 0.3,:) = NaN;
[R2,RMSE,NSE,PBIAS] = estimate_evaluation_metric(nanmean(gladseason,2)./100.*area,nanmean(swf_sea_cal,2).*area);
KGE = estimateKGE(nanmean(gladseason,2)./100.*area,nanmean(swf_sea_cal,2).*area);

str = {['\rho = ' num2str(round(sqrt(R2),2))] , ...
       ['NSE = ' num2str(round(NSE,2))]};
t = add_title(axs(1),str,18,'in');
t.Color = 'b';
t.Position(2) = t.Position(2) - 0.25;
set(axs(1),'Color',[0.8 0.8 0.8]);

axs(2) = axes('Position',[pos1(1) + pos1(3) + 0.02 pos1(2) 0.20 pos1(4)]);

sim = nanmean(swf_sea_cal,2).*area.*re^2./1e6;
obs = nanmean(gladseason,2)./100.*area.*re^2./1e6;
obs(isnan(sim)) = NaN;
sim2d = NaN(720,360);
obs2d = NaN(720,360);
sim2d(index_lnd) = sim;
obs2d(index_lnd) = obs;

plot(axs(2),nansum(obs2d,1),-89.75:0.5:89.75,'k-','LineWidth',1.5);hold on; grid on;
plot(axs(2),nansum(sim2d,1),-89.75:0.5:89.75,'r--','LineWidth',1.5); 
ylim([-60 80]);

leg = legend('GLAD','Simulation','FontSize',18,'FontWeight','bold');

[R2,RMSE,NSE,PBIAS,MSE,NSE1] = estimate_evaluation_metric(nansum(obs2d,1)',nansum(sim2d,1)');

str = {['\rho = ' num2str(round(sqrt(R2),2))] , ...
       ['NSE = ' num2str(round(NSE,2))]};
t = add_title(axs(2),str,18,'in');
t.Color = 'b';
t.Position(1) = t.Position(1) + 0.10;
t.Position(2) = t.Position(2) - 0.25;
add_title(axs(2),'(b) Latitude average',20,'out');

swf_yr_cal(lakein,:) = NaN; 
swf_yr_cal = swf_yr_cal(:,(2001 - 1993 + 1) : end);
continent = struct([]);
continent_code = {'af',    'ar',    'as',  'au',        'eu',    'gr',       'na',           'sa',           'si'     };
continent_name = {'Africa','Arctic','Asia','Austrialia','Europe','Greenland','North America','South America','Siberia'};

cols = distinguishable_colors(9);

wtot = axs(2).Position(3) + axs(2).Position(1) - axs(1).Position(1);
margin = 0.05;
w = (wtot - margin*3)/4;

k = 1;
%figure;

pw(frac < 1) = NaN;
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
numc = 0;
swf_yr_cal(isnan(gladannual)) = NaN;
swf_yr_cal  = swf_yr_cal  .* area;
gladannual  = gladannual./100 .* area;
gladseason  = gladseason./100 .* area;

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
   numc = numc + length(continent(i).index);
   %plot(xc(continent(i).index),yc(continent(i).index),'.','Color',cols(k,:)); hold on;

   if i <= 4
       axs(k+2) = axes('Position',[axs(1).Position(1) + (k-1)*(w+margin) axs(1).Position(2) - w - 0.1 w w]);
   else
       axs(k+2) = axes('Position',[axs(1).Position(1) + (k-5)*(w+margin) axs(1).Position(2) - 2*w - 0.15 w w]);
   end
   plot(axs(k+2),2001:2014,zscore(nansum(gladannual(continent(i).index,:),1)),'k-','LineWidth',2); hold on;grid on;
   plot(axs(k+2),2001:2014,zscore(nansum(swf_yr_cal(continent(i).index,:),1)),'r--','LineWidth',2);
   axs(k+2).GridLineStyle = '-';
   if i <= 4
       xlim([2001 2014]);
       set(gca,'XTick',[2002:2:2014])
       xticklabels("")
   else
       xlim([2001 2014]);
   end
   add_title(axs(k+2),[labels{k+1} ' ' continent(i).name],15,'out');
   k = k + 1;
end

ann = [];
sea = [];
sim = [];
for i = [1 2 3 4 5 7 8 9]
    ann = [ann; nansum(nanmean(gladannual(continent(i).index,:),2))];
    sea = [sea; nansum(nanmean(gladseason(continent(i).index,:),2))];
    sim = [sim; nansum(nanmean(swf_yr_cal(continent(i).index,:),2))];
end
figure;
%plot(sea,ann,'bx','LineWidth',2); hold on; grid on;
%plot(sim,ann,'ro','LineWidth',2);
plot(sim,sea,'ro','LineWidth',2);hold on; grid on;
plot([0 3e5],[0 3e5],'k-','LineWidth',2);



