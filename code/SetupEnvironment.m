function [e3sm_input, exportfig] = SetupEnvironment()
    
    home = getenv('HOME');
    exportfig = false;
    if ~isempty(getenv('exportfig')) && strcmp(getenv('exportfig'),'true')
        exportfig = true;
    end

    if strcmp(home,'/Users/xudo627')
        disp('Working on LOCAL!');

        addpath('/Users/xudo627/Developments/mylib/m/');
        addpath('/Users/xudo627/Developments/mylib/data/');
        addpath('/Users/xudo627/Developments/getPanoply_cMap/');
        % https://github.com/wschwanghart/topotoolbox
        addpath('/Users/xudo627/Developments/ofm_petsc/topotoolbox/');
        addpath('/Users/xudo627/Developments/ofm_petsc/topotoolbox/utilities/');
        addpath('/Users/xudo627/Developments/ofm_petsc/jigsaw-matlab/');
        addpath('/Users/xudo627/Developments/ofm_petsc/Matlab_Scripts/');
        % https://github.com/dengwirda/inpoly
        addpath('/Users/xudo627/Developments/inpoly/');
        % https://github.com/donghuix/Setup-E3SM-Mac
        addpath('/Users/xudo627/Developments/Setup-E3SM-Mac/matlab-scripts-to-process-inputs/');
        addpath('/Users/xudo627/Developments/Setup-E3SM-Mac/matlab-scripts-for-mosart/');
        % https://github.com/g2e/m_map
        addpath('/Users/xudo627/Developments/m_map/');
        addpath('/Users/xudo627/Developments/petsc/share/petsc/matlab/');
        
        e3sm_input = '/Users/xudo627/Library/CloudStorage/OneDrive-PNNL/projects/cesm-inputdata/';

    elseif strcmp(home,'/global/homes/d/donghui')
        disp('Working on NERSC!');
        addpath([home '/mylib/m/']);
        addpath([home '/mylib/data/']);
        addpath([home '/getPanoply_cMap/']);
        addpath([home '/Setup-E3SM-Mac/matlab-scripts-to-process-inputs/']);
        addpath([home '/Setup-E3SM-Mac/matlab-scripts-for-mosart/']);
        addpath([home '/inpoly/']);
        addpath([home '/m_map/']);

        e3sm_input = '/global/cfs/cdirs/e3sm/inputdata/';
        
    elseif strcmp(home,'/qfs/people/xudo627')
        disp('Working on COMPY!');
        
        addpath('/qfs/people/xudo627/mylib/m/');
        addpath('/qfs/people/xudo627/mylib/data/');
        addpath('/qfs/people/xudo627/getPanoply_cMap/');
        addpath('/qfs/people/xudo627/Setup-E3SM-Mac/matlab-scripts-to-process-inputs/');
        addpath('/qfs/people/xudo627/Setup-E3SM-Mac/matlab-scripts-for-mosart/');
        addpath('/qfs/people/xudo627/inpoly/');
        addpath('/qfs/people/xudo627/m_map/');
        
        exportfig = false; % exportgraphics is not available on compy matlab
        
        e3sm_input = '/compyfs/inputdata/';

    else
        disp('Need to install necessary packages!')
        stop(['Cannot work on ' home '!']);

    end

end
