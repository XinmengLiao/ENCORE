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
dlkcat_scriptDir = '/Users/xinmengliao/Documents/Project/20231127_GEMs/dlkcat';
minus_growth = table([], [], [], 'VariableNames', {'sample', 'bin', 'growth'});

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

% STEP 2 Load or import GEM models
    %model = importModel(model_file);
    model = loadConventionalGEM(); % Function has been changed locally, PDC haven't yet
    % check the reactions in ConventionalGEM
    [model.rxns,constructEquations(model)]
    %eqn=constructEquations(model)

% STEP 3 Convert the traditional one to ecGEMs one, but this ecGEMs is not functioning now;  
    % false is a full model, which is used in our pipeline 
    [ecModel] = makeEcModel(model,false,[],GeneSeqData); 

    % Backup the models
    %saveEcModel(ecModel,'ecGEMtest_stage1.yml');
    %model_tmp = fullfile(currentPath, bin, 'models','ecGEM_stage1.yml')
    %writeYAMLmodel(ecModel,  model_tmp)


%% STAGE 2: Integration of Kcat into the ecGEMs
% STEP 4 Gather EC numbers from Uniprot and KEGG, required for BRENDA as kcat source
    % Importmodel function in the loadEcModel() need to be set as false 
    %ecModel=loadEcModel('ecGEM_stage1.yml'); 
    
    % Uniprot EC number searching 
    kegg.ID      = params.kegg.ID;
    uniprot.ID   = params.uniprot.ID;
    filePath    = fullfile(params.path,'data');
    uniprotPath = fullfile(filePath,'uniprot.tsv');
    uniprot.type = params.uniprot.type;
    if strcmp(uniprot.type,'taxonomy')
        uniprot.type = 'taxonomy_id';
    end
    if params.uniprot.reviewed
        uniprotRev = 'reviewed:true+AND+';
    else
        uniprotRev = '';
    end

    
    % Download the uniprot database
    if ~strcmp(params.uniprot.ID, '')
        uniproturl = ['https://rest.uniprot.org/uniprotkb/stream?query=' uniprotRev ...
                   uniprot.type ':' num2str(uniprot.ID) '&fields=accession%2Cgene_names%2Cec%2Cmass%2Csequence&format=tsv&compressed=false&sort=protein_name%20asc'];
    
        try
            urlwrite(uniproturl,uniprotPath,'Timeout',30);
            fid         = fopen(uniprotPath,'r');
            fileContent = textscan(fid,'%q %q %q %q %q','Delimiter','\t','HeaderLines',1);
            if size(fileContent{1},1) > 1
                fclose(fid);
                fprintf('Model-specific Uniprot database stored at %s\n',uniprotPath);
                uniprot_exist = true;
                uniprot_notfound = false; 
            else
                fclose(fid);
                fprintf('Empty model-specific Uniprot database stored at %s\n',uniprotPath);
                uniprot_exist = false;
                uniprot_notfound = true; 
            end
        catch
            error(['Download failed, check your internet connection and try again, or manually download: ' uniproturl])
        end
    else
        uniprot_exist = false;
        uniprot_notfound = true; 
    end
    

    % Download data from KEGG
    if ~strcmp(params.kegg.ID, '')
        kegg_url = ['http://rest.kegg.jp/list/' kegg.ID];
        keggPath = fullfile(filePath,'kegg.tsv');
        try
            options = weboptions('Timeout', 30); 
            websave(keggPath, kegg_url, options);
            keggfid         = fopen(keggPath,'r');
            keggfileContent = textscan(keggfid,'%q %q %q %q %q %q %q','Delimiter',',','HeaderLines',0);
            if size(keggfileContent{1},1) > 1
                fclose(keggfid);
                fprintf('Model-specific KEGG database stored at %s\n',keggPath);
                kegg_exist = true;
                kegg_notfound = false; 
            else
                fclose(keggfid);
                fprintf('Empty model-specific Uniprot database stored at %s\n',keggPath);
                kegg_exist = false;
                kegg_notfound = true; 
            end
        catch
            error(['Download failed, check your internet connection and try again, or manually download: ' kegg_url])
        end
    else
        kegg_exist = false;
        kegg_notfound = true; 
    end

    
    % merge Uniprot and/or KEGG EC number together into model
    %if ~strcmp(params.kegg.ID, '') || ~strcmp(params.uniprot.ID, '')
    if uniprot_exist && kegg_exist
        ecModel         = getECfromGEM(ecModel);
        noEC = cellfun(@isempty, ecModel.ec.eccodes);
        ecModel         = getECfromDatabase(ecModel,noEC); % Gather Kcat from Brenda into model based on Uniprot and KEGG
        kcatList_fuzzy  = fuzzyKcatMatching(ecModel); % The fuzzyKcatMatching function has been changed locally, but the one in PDC still remaines unchanged. 
        kcatlist = true
    elseif uniprot_exist && ~kegg_exist
        ecModel         = getECfromGEM(ecModel);
        noEC = cellfun(@isempty, ecModel.ec.eccodes);
        ecModel         = getECfromDatabase_uniprot(ecModel,noEC); % Gather Kcat from Brenda into model based on Uniprot and KEGG
        kcatList_fuzzy  = fuzzyKcatMatching(ecModel); % The fuzzyKcatMatching function has been changed locally, but the one in PDC still remaines unchanged. 
        kcatlist = true
    elseif ~uniprot_exist && kegg_exist
        ecModel         = getECfromGEM(ecModel);
        noEC = cellfun(@isempty, ecModel.ec.eccodes);
        ecModel         = getECfromDatabase_kegg(ecModel,noEC); % Gather Kcat from Brenda into model based on Uniprot and KEGG
        kcatList_fuzzy  = fuzzyKcatMatching(ecModel); % The fuzzyKcatMatching function has been changed locally, but the one in PDC still remaines unchanged. 
        kcatlist = true
    else
        kcatlist = false

    end


% STEP 5 Gather metabolite SMILES, required for DLKcat as kcat source
    % Metabolite SMILES are gathered from PubChem.
    % Finding the chemical structures from PubChem for every metabolite
    [ecModel, noSmiles] = findMetSmiles(ecModel);

    %backup the model 
    %model_tmp = fullfile(currentPath, bin, 'models','ecGEM_stage2.yml')
    %writeYAMLmodel(ecModel,  model_tmp)
    %ecModel=loadEcModel('ecGEM_stage2.yml'); 

% STEP 6 Predict kcat values with DLKcat
    % Running DLKcat in Docker (very slow) or locally by python3.12 (miniconda3)
    writeDLKcatInput(ecModel,[],[],[],[],true);
    % runDLKcat();

    setenv('PATH', strcat('/Users/xinmengliao/miniconda3/bin:', getenv('PATH')));
    dlkcat_result_nokcat = fullfile(ecGEMpath, bin, 'data','DLKcat.tsv');
    dlkcat_results_kcat = fullfile(ecGEMpath, bin, 'data','DLKcatOutput.tsv');
    dlkcat_py = fullfile(dlkcat_scriptDir, 'DLKcat.py');
    cd(dlkcat_scriptDir)
    command = sprintf('python3 "%s" "%s" "%s"', dlkcat_py, dlkcat_result_nokcat, dlkcat_results_kcat);
    status = system(command);
    if status == 0 && exist(fullfile(params.path,'data/DLKcatOutput.tsv'))
        delete(fullfile(params.path,'/data/DLKcat.tsv'));
        movefile(fullfile(params.path,'/data/DLKcatOutput.tsv'), fullfile(params.path,'/data/DLKcat.tsv'));
        disp('DLKcat prediction completed.');
    else    
        error('DLKcat encountered an error or it did not create any output file.')
    end

    cd(fullfile(ecGEMpath,bin))

% STEP 7 Load DLKcat output
    kcatList_DLKcat = readDLKcatOutput(ecModel);

% STEP 8 Combine kcat from BRENDA and DLKcat
    if kcatlist
        kcatList_merged = mergeDLKcatAndFuzzyKcats(kcatList_DLKcat, kcatList_fuzzy)
    else
        kcatList_merged = kcatList_DLKcat;
    end
     

% STEP 9 Take kcatList and populate edModel.ec.kcat
    ecModel  = selectKcatValue(ecModel, kcatList_merged);

% STEP 10 Get kcat values across isozymes
    ecModel = getKcatAcrossIsozymes(ecModel);

% STEP 11 Calculate molecular weight and store them as model.ec.mw 
    ecModel.ec.mw = zeros(size(ecModel.ec.sequence));
    for i = 1:length(ecModel.ec.mw)
        if strcmp(ecModel.ec.sequence{i},'Sequence unavailable')
            ecModel.ec.mw(i) = 0;
        elseif isnan(ecModel.ec.sequence{i})
            ecModel.ec.mw(i) = 0;
        else
            try
                ecModel.ec.mw(i) = molweight(strrep(ecModel.ec.sequence{i},'X','D'));
            catch 
                fprintf('Protein No. %d does not have molecular weight. Mean molecular weight will be used.\n', i);
                ecModel.ec.mw(i) = 0;
            end
        end
    end
    % Replace the missing protein MWs with avg enzyme weight in this species
    ecModel.ec.mw(ecModel.ec.mw==0) = mean(ecModel.ec.mw>0);


% STEP 12 Get standard kcat (need to load Uniprot database)
    if uniprot_notfound == false
        [ecModel, rxnsMissingGPR, standardMW, standardKcat] = getStandardKcat(ecModel);
    else
        disp("Will not get standard kcat from Uniprot since proteome can not be found in Uniprot")
    end

% STEP 13 Apply kcat constraints from ecModel.ec.kcat to ecModel.S
    ecModel = applyKcatConstraints(ecModel);

% STEP 14 Set upper bound of protein pool
    Ptot  = ModelAdapter.params.Ptot;
    f     = ModelAdapter.params.f;
    sigma = ModelAdapter.params.sigma;
    ecModel = setProtPoolSize(ecModel,Ptot,f,sigma);

% STEP 15 Save EC Model
    finalid = strcat(binID, 'Ready_ecGEM.yml')
    model_tmp = fullfile(ecGEMpath, bin, 'models',finalid)
    writeYAMLmodel(ecModel,  model_tmp)

% STEP 16 Run FVA for exchange reactions only 
    % Add gurobi solver into the working path
    addpath(genpath('/Library/gurobi1201/macos_universal2/matlab'));

    % set RAVEN solver as gurobi
    setRavenSolver('gurobi')

    % Double confirm if RAVEN solver is gurobi
    checkInstallation

    % load ecGEM
    ecmodel_path = strcat(model_tmp)
    ecModel = readYAMLmodel(ecmodel_path)

    % check the reactions in ecModel ensure exchange reations are completed 
    [ecModel.rxns, constructEquations(ecModel)];
    
    % Extract the exchange reactions
    exchangeRxns = getExchangeRxns(ecModel);
    prot_idx = startsWith(exchangeRxns, 'prot');
    sink_idx = startsWith(exchangeRxns, 'sink');
    remove_idx = prot_idx | sink_idx;
    exchangeRxns = exchangeRxns(~remove_idx);
    
    % Create matrix to store LB and UB
    n_exc = length(exchangeRxns);
    LB_exc = zeros(n_exc, 1);
    UB_exc = zeros(n_exc, 1); 
    equations= cell(n_exc, 1);
    
    % FVA
    sensitivity = 0.95;
    %ecModel = setParam(ecModel,'ub',ecModel.rxns(ecModel.c>0),1);
    sol_opt = solveLP(ecModel);
    sol0 = solveLP(ecModel);
    objFlux = sol_opt.x(ecModel.c>0)

    if objFlux <= 0
        parts = strsplit(ecGEMpath, '/'); 
        sample = parts{end};  
        sample = string(sample);
        bin = string(bin);
        tmp = table(sample, bin, objFlux, 'VariableNames', {'sample', 'bin', 'growth'});            
        minus_growth = [minus_growth; tmp]
    else
        % run FVA for (sensitivity rate) of the maximun growth rate 
        bioModel = setParam(ecModel,'lb',ecModel.rxns(ecModel.c>0),sensitivity * objFlux);
        bioModel = setParam(bioModel,'ub',ecModel.rxns(ecModel.c>0),objFlux);
        [found, idx] = ismember(exchangeRxns, bioModel.rxns);
        
        for i = 1:length(exchangeRxns)
            i;
            iRxn = bioModel.rxns(idx(i));
            iModel = setParam(bioModel,'obj',iRxn,1);
            iSol = solveLP(iModel);
            
            try
                UB_exc(i) = iSol.x(ismember(iModel.rxns, iRxn));
            catch
                UB_exc(i) = sol0.x(logical(iModel.c));
                tmp = table(sample, bin, "No UB_exc", 'VariableNames', {'sample', 'bin', 'growth'}); 
                minus_growth = [minus_growth; tmp];
            end
            
            iModel = setParam(bioModel,'obj',iRxn,-1);
            iSol = solveLP(iModel);
            try
                LB_exc(i) = iSol.x(ismember(iModel.rxns, iRxn));
            catch
                LB_exc(i) = sol0.x(logical(iModel.c));
                tmp = table(sample, bin, "No LB_exc", 'VariableNames', {'sample', 'bin', 'growth'}); 
                minus_growth = [minus_growth; tmp];
            end
           
            equations{i} = constructEquations(bioModel, iRxn);
        end

        exchange_results = table(exchangeRxns, equations, LB_exc, UB_exc);
        fva_filename = strcat(binID, 'FVA.txt')
        outputfile = fullfile(ecGEMpath, bin, 'output',finalid)
        writetable(exchange_results, outputfile,'Delimiter', '\t', 'QuoteStrings', false, 'FileType', 'text');
    end  
end

minus_file = fullfile(ecGEMpath, 'minus_growth.txt')
writetable(minus_growth, minus_file,'Delimiter', '\t', 'QuoteStrings', false, 'FileType', 'text');

