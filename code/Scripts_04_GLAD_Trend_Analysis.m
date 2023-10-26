clear;close all;clc;

[e3sm_input, exportfig] = SetupEnvironment();

domain_file = [e3sm_input 'share/domains/domain.lnd.r05_oEC60to30v3.190418.nc'];
frac = ncread(domain_file,'frac'); frac = flipud(frac');
frac = flipud(frac');

inputdata = '/Users/xudo627/Library/CloudStorage/OneDrive-PNNL/projects/cesm-inputdata/';
dx = 0.5; dy = 0.5;
lon = -180 + dx/2 : dx : 180 - dx/2;
lat = 90 - dy/2 : -dy : -90 + dy/2;
[lon,lat] = meshgrid(lon,lat);

[~,~,area] = xc2xv(lon,lat,dx,dy,true);

tmp1 = NaN(20,1);
tmp2 = NaN(20,1);
tmp3 = NaN(20,1);
tmp4 = NaN(20,360);
tmp5 = NaN(360,720,20);
slopes = NaN(360,1);
slopes05 = NaN(360,720);
k = 1;
frac = ncread([inputdata 'share/domains/domain.lnd.r05_oEC60to30v3.190418.nc'],'frac');
frac = flipud(frac');

continent = struct([]);
continent_code = {'af',    'ar',    'as',  'au',        'eu',    'gr',       'na',           'sa',           'si'     };
continent_name = {'Africa','Arctic','Asia','Austrialia','Europe','Greenland','North America','South America','Siberia'};
load('colorblind_colormap.mat');
figure;
for i = 1 : length(continent_code)
   code = continent_code{i};
   continent(i).code  = code;
   continent(i).name  = continent_name{i};
   continent(i).index = []; 
   continent(i).swf   = NaN(20,1);
   S = shaperead(['../data/HydroBASINS/hybas_' code '_lev01-06_v1c/hybas_' code '_lev01_v1c.shp']);
   for j = 1 : length(S)
       tmp = inpoly2([lon(:) lat(:)],[S(j).X' S(j).Y']);
       tmp = find(tmp == 1);
       %plot(S(j).X,S(j).Y,'-','LineWidth',2);
       continent(i).index = [continent(i).index; tmp];
   end
   plot(lon(continent(i).index),lat(continent(i).index),'.','Color',colorblind(i,:),'LineWidth',3); hold on;
end
legend(continent_name,'FontSize',13,'FontWeight','bold');


for yr = 2001 : 2020
    disp(yr);
    load(['../data/GLAD05/GLAD05_' num2str(yr) '.mat']);
    GLAD05(frac < 1) = NaN;
    
    ind1 = find(lon >= +55 & lon <= +60 & lat >= +50 & lat <= +55);
    ind2 = find(lon >= -80 & lon <= -60 & lat >= +50 & lat <= +60);
    ind3 = find(lat >= 40);
    tmp1(k) = nansum(GLAD05(ind1) .*area(ind1)) ./ nansum(area(ind1));
    tmp2(k) = nansum(GLAD05(ind2) .*area(ind2)) ./ nansum(area(ind2));
    tmp3(k) = nansum(GLAD05(ind3) .*area(ind3)) ./ nansum(area(ind3));
    tmp4(k,:) = nanmean(GLAD05./100.*area,2)./1e6;
    tmp5(:,:,k) = GLAD05;

    for ii = 1 : length(continent)
        tmp = nansum(GLAD05(continent(ii).index).*area(continent(ii).index)) ./ 100 ./ 1e6;
        continent(ii).swf(k) = tmp;
    end
    k = k + 1;
end
for i = 1 : 360
    disp(i);
    [slopes(i,1)] = sens_slope(tmp4(:,i));
    [z, trend, h, p(i,1)] = mk_test(tmp4(:,i));
end

k = 1;

for i = 2001 : 2020
    for j = 1 : 12
        yrs(k) = i;
        mos(k) = j;
        k = k + 1;
    end
end

for i = 1 : 360
    disp(i);
    for j = 1 : 720
        tmp = tmp5(i,j,:);
        tmp = tmp(:);
        if ~all(isnan(tmp))
            [slopes05(i,j)] = sens_slope(tmp);
        end
    end
end

figure;
subplot(1,3,1);
plot(tmp1,'r-','LineWidth',2); grid on;
subplot(1,3,2);
plot(tmp2,'g-','LineWidth',2); grid on;
subplot(1,3,3);
plot(tmp3,'k-','LineWidth',2); grid on;

figure;
[cb,lon,lat] = plot_globalspatial(lon,lat,slopes05,1,1);
clim([-0.25 0.25]); colormap(blue2red(121));
ylim([-60 80]);
pos = get(gca,'Position');
cb.Position(1) = pos(1)+pos(3) + 0.01;
cb.Position(2) = pos(2);
cb.Position(3) = 0.02;
cb.Position(4) = pos(4);
cb.AxisLocation = 'out';

figure;
for i = 1 : 9
    subplot(3,3,i);
    plot(continent(i).swf,'-','Color',colorblind(i,:),'LineWidth',2); hold on; grid on;
    add_title(gca,continent(i).name,15,'in');
end

T1 = ncread('~/DATA/CRU/cru_ts4.07.2001.2010.tmp.dat.nc','tmp');
T2 = ncread('~/DATA/CRU/cru_ts4.07.2011.2020.tmp.dat.nc','tmp');
P1 = ncread('~/DATA/CRU/cru_ts4.07.2001.2010.pre.dat.nc','pre');
P2 = ncread('~/DATA/CRU/cru_ts4.07.2011.2020.pre.dat.nc','pre');
Tmo  = cat(3,T1,T2); clear T1 T2;
Pmo  = cat(3,P1,P2); clear P1 P2;
Tyr  = NaN(360,720,20);
Pyr  = NaN(360,720,20);
k = 1;
tmp6 = NaN(20,360);
tmp7 = NaN(20,360);
for i = 2001 : 2020
    tmp = nanmean(Tmo(:,:,yrs == i),3);
    Tyr(:,:,k) = flipud(tmp');
    tmp = flipud(tmp');
    tmp(frac < 1) = NaN;
    tmp6(k,:) = nanmean(tmp,2);

    tmp = nanmean(Pmo(:,:,yrs == i),3);
    tmp = flipud(tmp');
    tmp(frac < 1) = NaN;
    tmp7(k,:) = nanmean(tmp,2);
    k = k + 1;
end

for i = 1 : 360
    disp(i);
    [slopesT(i,1)] = sens_slope(tmp6(:,i));
    [slopesP(i,1)] = sens_slope(tmp7(:,i));
    [z, trend, h, pT(i,1)] = mk_test(tmp6(:,i));
    [z, trend, h, pP(i,1)] = mk_test(tmp7(:,i));
end

figure;
dy = 0.5;
subplot(1,3,1);
lats = 90-dy/2 : -dy : -90 + dy/2;
plot(slopes,lats,'kx','LineWidth',2); hold on; grid on;
plot(slopes(p <= 0.05),lats(p <= 0.05),'rx','LineWidth',2); 
ind = find(~isnan(slopes));
p1 = polyfit(lats(ind),slopes(ind),2);
y1 = polyval(p1,lats(ind));
plot(y1,lats(ind),'b-','LineWidth',5); hold on; grid on;

subplot(1,3,2);
plot(slopesT,lats,'kx','LineWidth',2); hold on; grid on;
plot(slopesT(pT <= 0.05),lats(pT <= 0.05),'rx','LineWidth',2); 
ylim([-60 80]);
subplot(1,3,3);
plot(slopesP,lats,'kx','LineWidth',2); hold on; grid on;
plot(slopesP(pP <= 0.05),lats(pP <= 0.05),'rx','LineWidth',2); 
ylim([-60 80]);
