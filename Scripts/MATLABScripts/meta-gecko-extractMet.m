% Add GECKO and RAVEN path into MatLab path
%addpath(genpath('/cfs/klemming/projects/snic/naiss2023-23-637/ecgem/GECKO'));
%addpath(genpath('/cfs/klemming/projects/snic/naiss2023-23-637/ecgem/RAVEN'));

% Check installation
%GECKOInstaller.install
% create a GECKO project, with the default folders setting
%startGECKOproject(filename, currentPath)
% will create the folders in current path,  need to go inside and change
% the Adapter file 

%% Stage 1 Expansion from the traditional GEMs to ecGEMs
% STEP 1 Set modelAdapter path # store all files in ecModel_files 

ecGEMpath = ''
GEMpath = '';
cd(GEMpath)
GEMfiles = dir(fullfile(GEMpath, '*.xml'));
currentPath = pwd

for i = 1:length(GEMfiles)
    [~, bin, ~] = fileparts(GEMfiles(i).name);  % get bins name
    disp(['Current bin: ', bin]);
    binID = regexp(bin, "\.","split");
    binID = strjoin(binID(1:3),"" );
    adapter_file = strcat(binID, "ecGECKOAdapter.m");
    adapterLocation = strcat(ecGEMpath, "/",bin,"/", adapter_file);
    addpath(genpath(fullfile(ecGEMpath,bin)));
    ModelAdapter = ModelAdapterManager.setDefault(adapterLocation);
    ModelAdapter = ModelAdapterManager.getDefault();
    params = ModelAdapter.getParameters();
    cd(fullfile(ecGEMpath, bin));
    binpath= pwd ;

% load the protein sequences
    faa = dir(fullfile(binpath, '*.faa.txt'));
    faa = string({faa.name});
    GeneSeqData_file = fullfile(binpath, faa)
    GeneSeqData = readtable(GeneSeqData_file, 'Delimiter', '\t', 'ReadVariableNames', false);
    GeneSeqData.Properties.VariableNames = {'Gene', 'Protein_seq'};

    model = loadConventionalGEM(); 
    [model.rxns,constructEquations(model)]

    [ecModel] = makeEcModel(model,false,[],GeneSeqData); 
    
    [uniqueNames, ~, uniqueIdx] = unique(regexprep(ecModel.metNames,'^prot_.*',''));
    uniqueSmiles(1:numel(uniqueNames),1) = {''};
    metFile = strcat(binID,'.unique_metabolites.txt');
    metPath = fullfile(ecGEMpath, bin,"data",metFile)
    fid = fopen(metPath, 'w');

    for i = 1:length(uniqueNames)
        fprintf(fid, '%s\n', uniqueNames{i});
    end
    fclose(fid);

end
