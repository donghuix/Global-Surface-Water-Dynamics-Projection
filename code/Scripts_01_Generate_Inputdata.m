clear;close all;clc;

tag  = 'GLOBE';
ntot = 10;

[e3sm_input, exportfig] = SetupEnvironment();
re = 6.37122e6;% Earth radius

domain_file = [e3sm_input 'share/domains/domain.lnd.r05_oEC60to30v3.190418.nc'];
surf_file   = [e3sm_input 'lnd/clm2/surfdata_map/surfdata_0.5x0.5_simyr2010_c200624.nc'];
%mosart_file = [e3sm_input 'rof/mosart/MOSART_Global_half_20200720.nc'];
mosart_file = '../data/MOSART_Global_half_20200720.nc';

out_netcdf_dir = '../inputdata';

frac = ncread(domain_file,'frac');
mask = ncread(domain_file,'mask');
lonc = ncread(domain_file,'xc');
latc = ncread(domain_file,'yc');
area = ncread(domain_file,'area');

index_lnd = find(frac > 0 & latc >= -60);
lon1d = lonc(index_lnd);
lat1d = latc(index_lnd);
frac1d = frac(index_lnd);
mask1d = mask(index_lnd);
area1d = area(index_lnd);
[xv,yv] = xc2xv(lon1d,lat1d,0.5,0.5);

continent = struct([]);
continent_code = {'af',    'ar',    'as',  'au',        'eu',    'gr',       'na',           'sa',           'si'     };
continent_name = {'Africa','Arctic','Asia','Austrialia','Europe','Greenland','North America','South America','Siberia'};
cmap = getPanoply_cMap('EO_aura_omi_formal');
%cb = plot_globalspatial(lonc,latc,frac',1,1); colormap(cmap); hold on;
load('colorblind_colormap.mat');

%figure;
%for i = 1 : length(continent_code)
%    code = continent_code{i};
%    continent(i).code  = code;
%    continent(i).name  = continent_name{i};
%    continent(i).index = []; 
%    S = shaperead(['../data/HydroBASINS/hybas_' code '_lev01-06_v1c/hybas_' code '_lev01_v1c.shp']);
%    for j = 1 : length(S)
%        tmp = inpoly2([lon1d lat1d],[S(j).X' S(j).Y']);
%        tmp = find(tmp == 1);
        %plot(S(j).X,S(j).Y,'-','LineWidth',2);
%        continent(i).index = [continent(i).index; tmp];
%    end
%    plot(lon1d(continent(i).index),lat1d(continent(i).index),'.','Color',colorblind(i,:),'LineWidth',3); hold on;
%end
%legend(continent_name,'FontSize',13,'FontWeight','bold');

numc = length(lon1d);
% (1). Generate defualt MOSART input file in 1D format
usrdat_name = [tag '_1d'];
fname_out1 = CreateMOSARTUgrid(index_lnd, mosart_file, out_netcdf_dir, usrdat_name,1);

% figure;
% show_river_network(fname_out1,0.25,'b-');

%(2). Generate defualt domain file in 1D format
fname_out2 =  sprintf('%s/domain_lnd_%s%s.nc',out_netcdf_dir,tag,'_1d');
area = generate_lnd_domain(lon1d,lat1d,xv,yv,frac1d,mask1d,area1d,fname_out2);

%(3). Generate default ELM surface input file in 1D format
sand = NaN(numc,10);
clay = NaN(numc,10);
sandall = ncread(surf_file,'PCT_SAND');
clayall = ncread(surf_file,'PCT_CLAY');
for i = 1 : 10
    tmp = sandall(:,:,i);
    sand(:,i) = tmp(index_lnd);
    tmp = clayall(:,:,i);
    clay(:,i) = tmp(index_lnd);
end
clear sandall clayall;

fmax = ncread(surf_file,'FMAX');
fmax = fmax(index_lnd);
silt = 100 - sand - clay;
bsw_mu    = 2.91 + 0.159.*clay;
sucsat_mu = 10 .* ( 10.^(1.88-0.0131.*sand) );
xksat_mu  = 0.0070556 .*( 10.^(-0.884+0.0153.*sand) );
watsat_mu = 0.01.*( 48.9 - 0.126.*sand );

fname_out3 = CreateCLMUgridSurfdatForE3SM(...
                    index_lnd, surf_file,out_netcdf_dir, [tag '_t1d'],  ...
                    ones(numc,1).*2.5, ones(numc,1).*5.5e-3,            ... % fdrain,max_drain,
                    ones(numc,1).*6, ones(numc,1),                      ... % ice_imped,snoalb_factor,
                    ones(numc,1).*0.5, fmax,                            ... % fover,fmax
                    bsw_mu,sucsat_mu, xksat_mu,watsat_mu,               ... % bsw,sucsat,xksat,watsat,
                    ones(numc,1).*0.4, ones(numc,1).*0.14,              ... % fc, mu
                    [], []);           

%(4). Generate calibration MOSART file with 10 global 1d domains
areaTotal2 = ncread(fname_out1,'areaTotal2');
rdep_uq  = NaN(numc*ntot,1);
rwid_uq  = NaN(numc*ntot,1);

fname_out4 = sprintf('%s/MOSART_%s.nc',out_netcdf_dir,[tag '_cal']);
create_mosart_cal(fname_out1,fname_out4,numc,ntot)

%(5). Generate calibration domain file with 10 global 1d domains
fname_out5 = sprintf('%s/domain_lnd_%s.nc',out_netcdf_dir,[tag '_cal']);
generate_lnd_domain_uq(fname_out2,fname_out5,ntot);

%(6). Generate calibration surface dataset with 10 global 1d domains
rng(4); % Random seed
X = lhsdesign(100,2,'Criterion','correlation'); % fover, fc
fover = exp(log(0.1) + (log(5)-log(0.1))*X(:,1));
%fover = 0.1 + (0.5-0.1)*X(:,1);
fc    = 0.001 + (0.4-0.001)*X(:,2);
figure;
plot(fc,fover,'bx','LineWidth',2); hold on;

fdrai_uq       = ones(numc*ntot,1).*2.5;
fdqraimax_uq   = ones(numc*ntot,1).*5.5e-3;
ice_imped_uq   = ones(numc*ntot,1).*6;
fover_uq       = ones(numc*ntot,1).*0.5;
fmax_uq        = repmat(fmax,ntot,1);
bsw_uq         = repmat(bsw_mu,ntot,1);
sucsat_uq      = repmat(sucsat_mu,ntot,1);
xksat_uq       = repmat(xksat_mu,ntot,1);
watsat_uq      = repmat(watsat_mu,ntot,1);
fc_uq          = ones(numc*ntot,1).*0.4;
mu_uq          = ones(numc*ntot,1).*0.14;
micro_sigma_uq = NaN(numc*ntot,1);
kh2osfc_uq     = NaN(numc*ntot,1);

for i = 1 : 10
    if i >= 10
        fname_out6 = sprintf('%s/surfdata_%s.nc',out_netcdf_dir,[tag '_cal_' num2str(i)]);
    else
        fname_out6 = sprintf('%s/surfdata_%s.nc',out_netcdf_dir,[tag '_cal_0' num2str(i)]);
    end
    
    disp([' Generating ' fname_out6]);
    
    for j = 1 : ntot
        fover_uq((j-1)*numc+1:j*numc) = fover((i-1)*ntot+j);
        fc_uq((j-1)*numc+1:j*numc)    = fc((i-1)*ntot+j);
    end

    generate_lnd_surface_uq(fname_out3,fname_out6,ntot,          ...
                            fdrai_uq,fdqraimax_uq,               ... % fdrain,max_drain,
                            ice_imped_uq,ones(numc*ntot,1),      ... % ice_imped,snoalb_factor,
                            fover_uq, fmax_uq,                   ... % fover,fmax
                            bsw_uq,sucsat_uq,xksat_uq,watsat_uq, ... % bsw,sucsat,xksat,watsat,
                            fc_uq,mu_uq,                         ... % fc, mu
                            micro_sigma_uq,kh2osfc_uq);
end
