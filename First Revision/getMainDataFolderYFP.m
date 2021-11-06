function mainDataFolder = getMainDataFolderYFP (varargin)
mainDataFolder = 'Y:\homes\navvab.afrashteh\CloudStation\Navvab-Sam\Data\NeuroPhotonic Paper\YFP Data';

if ~exist(mainDataFolder,'dir')
    mainDataFolder = '\\mohajerani-nas.uleth.ca\storage2\homes\navvab.afrashteh\CloudStation\Navvab-Sam\Data\NeuroPhotonic Paper\YFP Data';
end
return;

if nargin == 1
    mainDataFolder = 'E:\Users\samsoon.inayat\ProcessedData';
    return;
end
mainDataFolder = 'Y:\homes\navvab.afrashteh\CloudStation\Navvab-Sam\Data\NeuroPhotonic Paper\YFP Data';
if ~exist(mainDataFolder,'dir')
    mainDataFolder = 'F:\CloudStation\Navvab-Sam\Data\NeuroPhotonic Paper\YFP Data';
end
if ~exist(mainDataFolder,'dir')
    mainDataFolder = 'E:\Users\navvab.afrashteh\CloudStation\Navvab-Sam\Data\NeuroPhotonic Paper\YFP Data';
end
if ~exist(mainDataFolder,'dir')
    mainDataFolder = 'E:\Users\samsoon.inayat\CloudStation\Navvab-Sam\Data\NeuroPhotonic Paper\YFP Data';
end
