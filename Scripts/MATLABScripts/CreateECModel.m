% Check installation
addpath(genpath('/cfs/klemming/projects/snic/naiss2023-23-637/ecgem/GECKO'));
addpath(genpath('/cfs/klemming/projects/snic/naiss2023-23-637/ecgem/RAVEN'));
GECKOInstaller.install
% create a GECKO project, with the default folders setting
filename = ''
currentPath = ''
startGECKOproject(filename, currentPath)
