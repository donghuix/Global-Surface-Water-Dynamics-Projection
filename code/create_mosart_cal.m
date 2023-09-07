function create_mosart_cal(fin,fout,numc,ntot)

ncid_inp = netcdf.open(fin,'NC_NOWRITE');
ncid_out = netcdf.create(fout,'NC_CLOBBER');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid_inp);
dimid(1) = netcdf.defDim(ncid_out,'gridcell',ntot*numc);
dimid(2) = netcdf.defDim(ncid_out,'nele',11);

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
%                           Define variables
%
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
for ivar = 1 : nvars
    [varname,xtype,dimids,natts] = netcdf.inqVar(ncid_inp,ivar-1);
    switch varname
        case {'ele'}
            varid(ivar) = netcdf.defVar(ncid_out,varname,xtype,dimids); 
        otherwise
            varid(ivar) = netcdf.defVar(ncid_out,varname,xtype,dimid(1)); 
    end
    varnames{ivar} = varname;
    for iatt = 1 : natts
        attname = netcdf.inqAttName(ncid_inp,ivar-1,iatt-1);
        attvalue = netcdf.getAtt(ncid_inp,ivar-1,attname);
        
        netcdf.putAtt(ncid_out,ivar-1,attname,attvalue);
    end
end
varid = netcdf.getConstant('GLOBAL');

[~,user_name]=system('echo $USER');
netcdf.putAtt(ncid_out,varid,'Created_by' ,user_name(1:end-1));
netcdf.putAtt(ncid_out,varid,'Created_on' ,datestr(now,'ddd mmm dd HH:MM:SS yyyy '));
netcdf.endDef(ncid_out);

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
%                           Copy variables
%
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
for ivar = 1:nvars
    [varname,vartype,vardimids,varnatts]=netcdf.inqVar(ncid_inp,ivar-1);
    data = netcdf.getVar(ncid_inp,ivar-1);
    switch varname
        case {'ID','dnID'}
            data1 = data;
            for i = 2 : ntot
                tmp = data1;
                for j = 1 : length(tmp)
                    if tmp(j) ~= -9999
                        tmp(j) = tmp(j) + (i-1)*numc;
                    end
                end
                data = [data; tmp];
            end
        case {'rslp'}
            data(data < 0) = 0;
            data = repmat(data,ntot,1);
        case {'tslp'}
            data(data < 0) = 0.0001;
            data = repmat(data,ntot,1);
        case {'hslp'}
            data(data < 0) = 0.005;
            data = repmat(data,ntot,1);
        otherwise
            data = repmat(data,ntot,1);
    end
    netcdf.putVar(ncid_out,ivar-1,data);
end


netcdf.close(ncid_inp);
netcdf.close(ncid_out);

end