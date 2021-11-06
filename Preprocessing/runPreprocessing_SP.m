function runPreprocessing_SP (animalNumber,varargin)

Rotate = 0;
theta = 0;
Glut = 0;
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
    sFolders = dir(sprintf('%s\\*.tif',dataFolder));
    for jj = 1:length(sFolders)
        temp = sFolders(jj).name;
        if isempty(strfind(lower(temp),'spon'))
            continue;
        end
        PreProcessSpon(sprintf('%s\\',dataFolder),temp(1:(end-4)),Rotate,theta,Glut,AC);
    end
end
