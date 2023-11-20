clear;close all;clc;

[e3sm_input, exportfig] = SetupEnvironment();

yr = 2001;
load('../data/cal_9.mat','tsa');
xv = ncread('../data/domain_lnd_GLOBE_1d.nc','xv');
yv = ncread('../data/domain_lnd_GLOBE_1d.nc','yv');

load coastlines;

figure; set(gcf,'Position',[10 10 1200 600]);
patch(xv,yv,nanmean(tsa(1:66924,(2001-1993)*12+1:(2001-1993+1)*12)-273.15,2),'LineStyle','none'); hold on;
plot(coastlon,coastlat,'k-','LineWidth',1); grid on;
colormap(jet); colorbar; clim([-10 24]);
xlim([65 110]);
ylim([25 45]);

S = shaperead('../data/TPBoundary_new(2021)/TPBoundary_new(2021).shp');
plot(S.X,S.Y,'k-','LineWidth',2);

files = dir('~/Downloads/MOD11C3*.hdf');
lst = NaN(3600,7200,12);
for i = 1 : length(files)
    filename = fullfile(files(i).folder,files(i).name);
    disp(filename);
    lst_day    = double(hdfread(filename,'LST_Day_CMG'));
    lst_night  = double(hdfread(filename,'LST_Night_CMG'));
    day_time   = double(hdfread(filename,'Day_view_time'));
    night_time = double(hdfread(filename,'Night_view_time'));
    
    lst_day(lst_day == 0) = NaN;
    lst_day = lst_day.*0.02;
    lst_night(lst_night == 0) = NaN;
    lst_night = lst_night.*0.02;
    
    day_time(day_time == 255) = NaN;
    day_time = day_time .* 0.2;
    night_time(night_time == 255) = NaN;
    night_time = night_time .* 0.2;
    
    lst(:,:,i) = lst_day .* day_time./(day_time + night_time) + lst_night .* night_time./(day_time + night_time);
end

lon = -179.975 : 0.05 : 179.975;
lat = 89.975 : -0.05 : -89.975;
[lon,lat] = meshgrid(lon,lat);

figure; set(gcf,'Position',[10 10 1200 600]);
for i = 1 : 12
    subplot(3,4,i)
    imagesc([lon(1,1) lon(end,end)],[lat(1,1) lat(end,end)],lst(:,:,i)-273.15); hold on;
    set(gca,'YDir','normal'); 
    plot(coastlon,coastlat,'k-','LineWidth',1); grid on;
    colormap(jet); colorbar; clim([-10 24]);
    title(['MODIS: 2001-' num2str(i)],'FontSize',13,'FontWeight','bold');
    xlim([65 110]);
    ylim([25 45]);
end

figure; set(gcf,'Position',[10 10 1200 600]);
for i = 1 : 12
    tmp = lst(:,:,i)-273.15;
    tmp = convert_res(tmp,1,10)./100;
    subplot(3,4,i)
    imagesc([-179.75 179.75],[89.75 -89.75],tmp); hold on;
    set(gca,'YDir','normal'); 
    plot(coastlon,coastlat,'k-','LineWidth',1); grid on;
    colormap(jet); colorbar; clim([-10 24]);
    title(['Upscaled MODIS: 2001-' num2str(i)],'FontSize',13,'FontWeight','bold');
    xlim([65 110]);
    ylim([25 45]);
end


sim = tsa(1:66924,(2001-1993)*12+1:(2001-1993+1)*12)-273.15;
figure; set(gcf,'Position',[10 10 1200 600]);
for i = 1 : 12
    subplot(3,4,i)
    patch(xv,yv,sim(:,i),'LineStyle','none'); hold on;
    plot(coastlon,coastlat,'k-','LineWidth',1); grid on;
    colormap(jet); colorbar; clim([-10 24]);
    title(['E3SM: 2001-' num2str(i)],'FontSize',13,'FontWeight','bold');
    xlim([65 110]);
    ylim([25 45]);
end