function runPreprocessing_EV (animalNumber,varargin)
Rotate = 0;
theta = 0;
Glut = 0;
AC = 0;
if nargin > 1
    Rotate = varargin{1};
end
if Rotate && nargin > 2
    theta = varargin{2};
end
if nargin > 3    
    if strcmpi(varargin{3},'glut')
        mainDataFolder = getMainDataFolderGlut;
        Glut = 1;
    elseif strcmpi(varargin{3},'ac')
        mainDataFolder = getMainDataFolderAC;
        AC = 1;
    end
else
    mainDataFolder = getMainDataFolder;
end

for an = 1:length(animalNumber)
    dataFolder = makeName(animalNumber{an},mainDataFolder);
    eDataFolder = makeName('Evoked',dataFolder);
    folders = dir(eDataFolder);
    peDataFolder = makeName('pEvoked',dataFolder);
    for ii = 1:length(folders)
        if strcmp(folders(ii).name,'.') || strcmp(folders(ii).name,'..')
            continue;
        end
        filePath = makeName(folders(ii).name,eDataFolder);
        sFolders = dir(sprintf('%s\\*.tif',filePath));
        trial = 1;
        for jj = 1:length(sFolders)
            thisName = sFolders(jj).name;
            thisName = thisName(1:(end-4));
            if (strcmp(thisName(1:2),'hi') || strcmp(thisName(1:2),'lo') || strcmp(thisName(1:2),'t1'))
                if (strcmp(thisName(1:2),'t1') && ~AC)
                    continue;
                end
                if ~AC && ~Glut && ~exist([filePath,'\no',thisName(3:end),'.tif'])
                    continue;
                end
                file = thisName(1:end);
                savePath = makeName(folders(ii).name,peDataFolder);
                savePath = makeName(sprintf('trial_%.2d\\',trial),savePath);
                if ~exist(savePath,'dir')
                    mkdir(savePath);
                end
                PreProcessEV_HiLo(file,[filePath '\'],savePath,Rotate,theta,Glut,AC);
                trial = trial + 1;
            end
        end
    end
end
