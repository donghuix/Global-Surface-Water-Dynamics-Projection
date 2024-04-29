function [lakein,lakein2d] = getLakeIndex(e3sm_input)
    domain_file = [e3sm_input 'share/domains/domain.lnd.r05_oEC60to30v3.190418.nc'];

    load('LargeLakes.mat');
    xc   = ncread('../data/domain_lnd_GLOBE_1d.nc','xc');
    yc   = ncread('../data/domain_lnd_GLOBE_1d.nc','yc');
    lon  = ncread(domain_file,'xc');
    lat   = ncread(domain_file,'yc');
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

end