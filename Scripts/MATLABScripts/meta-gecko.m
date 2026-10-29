% Add GECKO and RAVEN path into MatLab path

%% Stage 1 Expansion from the traditional GEMs to ecGEMs
% STEP 1 Set modelAdapter path # store all files in ecModel_files 

ecGEMpath = ''
GEMpath = '';
cd(ecGEMpath);
% Only the bins with taxonomic information have a folder in ecGEMpath
GEMfiles = dir(fullfile(ecGEMpath, 'bin*'));
GEMfiles = GEMfiles([GEMfiles.isdir]);
currentPath = pwd;
parts = strsplit(ecGEMpath, '/');
sample = string(parts{end});
growth_table = table(strings(0,1), strings(0,1), zeros(0,1), 'VariableNames', {'sample', 'bin', 'growth'});

for i = 1:length(GEMfiles)
    bin = GEMfiles(i).name;  % get bins name
    disp(['Current bin: ', bin]);
    binID = regexp(bin, "\.","split");
    binID = strjoin(binID(1:3),"" );
    % Resume: a bin that finished in an earlier run (its growth file is written at the
    % very end of the bin) is skipped, and its growth is read back
    binGrowthFile = fullfile(ecGEMpath, bin, 'output', strcat(sample, ".", binID, 'growth.txt'));
    if exist(binGrowthFile, 'file')
        disp(strcat("Bin ", bin, " already finished in an earlier run, skipped."));
        fid = fopen(binGrowthFile, 'r');
        c = textscan(fid, '%s %s %f', 'Delimiter', '\t');
        fclose(fid);
        growth_table = [growth_table; table(string(c{1}), string(c{2}), c{3}, 'VariableNames', {'sample', 'bin', 'growth'})];
        continue
    end
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
    GeneSeqData_file = fullfile(binpath, faa);
    GeneSeqData = readtable(GeneSeqData_file, 'Delimiter', '\t', 'ReadVariableNames', false);
    GeneSeqData.Properties.VariableNames = {'Gene', 'Protein_seq'};

% STEP 2 Load or import GEM models
    model = loadConventionalGEM(); % replace this function by the one in Scripts/InstallationScripts/loadConventionalGEM.m before running. 
    [model.rxns,constructEquations(model)];

% STEP 3 Convert the traditional one to ecGEMs one, but this ecGEMs is not functioning now;  
    [ecModel] = makeEcModel(model,false,[],GeneSeqData); 


%% STAGE 2: Integration of Kcat into the ecGEMs
% STEP 4 Gather EC numbers from Uniprot and KEGG, required for BRENDA as kcat source

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
            for attempt = 1:6
                try
                    urlwrite(uniproturl,uniprotPath,'Timeout',30);
                    break
                catch downloadErr
                    if attempt == 6, rethrow(downloadErr); end
                    pause(2^attempt); % UniProt sometimes refuses bursts of requests: wait and retry
                end
            end
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
            error(['Download failed, check your internet connection and try again, or manually download: ' uniproturl]);
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
            for attempt = 1:6
                try
                    websave(keggPath, kegg_url, options);
                    break
                catch downloadErr
                    if attempt == 6, rethrow(downloadErr); end
                    pause(2^attempt); % KEGG sometimes refuses bursts of requests: wait and retry
                end
            end
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
            error(['Download failed, check your internet connection and try again, or manually download: ' kegg_url]);
        end
    else
        kegg_exist = false;
        kegg_notfound = true; 
    end

    
    % merge Uniprot and/or KEGG EC number together into model
    if uniprot_exist && kegg_exist
        ecModel         = getECfromGEM(ecModel);
        noEC = cellfun(@isempty, ecModel.ec.eccodes);
        ecModel         = getECfromDatabase(ecModel,noEC); % Gather Kcat from Brenda into model based on Uniprot and KEGG
        kcatList_fuzzy  = fuzzyKcatMatching(ecModel); % The fuzzyKcatMatching function has been changed locally, but the one in PDC still remaines unchanged. 
        kcatlist = true;
    elseif uniprot_exist && ~kegg_exist
        ecModel         = getECfromGEM(ecModel);
        noEC = cellfun(@isempty, ecModel.ec.eccodes);
        ecModel         = getECfromDatabase_uniprot(ecModel,noEC); % Gather Kcat from Brenda into model based on Uniprot and KEGG
        kcatList_fuzzy  = fuzzyKcatMatching(ecModel); % The fuzzyKcatMatching function has been changed locally, but the one in PDC still remaines unchanged. 
        kcatlist = true;
    elseif ~uniprot_exist && kegg_exist
        ecModel         = getECfromGEM(ecModel);
        noEC = cellfun(@isempty, ecModel.ec.eccodes);
        ecModel         = getECfromDatabase_kegg(ecModel,noEC); % Gather Kcat from Brenda into model based on Uniprot and KEGG
        kcatList_fuzzy  = fuzzyKcatMatching(ecModel); % The fuzzyKcatMatching function has been changed locally, but the one in PDC still remaines unchanged. 
        kcatlist = true;
    else
        kcatlist = false;

    end


% STEP 5 Gather metabolite SMILES, required for DLKcat as kcat source
    % Metabolite SMILES are gathered from PubChem.
    % Finding the chemical structures from PubChem for every metabolite
    % SMILES are only used to write the DLKcat input: skip if data/DLKcat.tsv already exists
    if ~exist(fullfile(params.path,'data/DLKcat.tsv'), 'file')
        [ecModel, noSmiles] = findMetSmiles(ecModel);
    else
        disp('DLKcat.tsv already exists, SMILES search skipped.');
    end


% STEP 6 Predict kcat values with DLKcat
    % Skip the prediction if data/DLKcat.tsv (with predicted kcats) already exists
    if ~exist(fullfile(params.path,'data/DLKcat.tsv'), 'file')
        % Running DLKcat in Docker (very slow) or locally by python3.12 (miniconda3)
        writeDLKcatInput(ecModel,[],[],[],[],true);
        try
            runDLKcat();
        catch dlkcatErr
            % Remove the input-only DLKcat.tsv (no kcat column), otherwise a rerun would skip prediction
            if exist(fullfile(params.path,'data/DLKcat.tsv'), 'file')
                delete(fullfile(params.path,'data/DLKcat.tsv'));
            end
            rethrow(dlkcatErr);
        end
        disp('DLKcat prediction completed.');
    else
        disp('DLKcat.tsv already exists, DLKcat prediction skipped.');
    end
    cd(fullfile(ecGEMpath,bin));

% STEP 7 Load DLKcat output
    kcatList_DLKcat = readDLKcatOutput(ecModel);

% STEP 8 Combine kcat from BRENDA and DLKcat
    if kcatlist
        kcatList_merged = mergeDLKcatAndFuzzyKcats(kcatList_DLKcat, kcatList_fuzzy);
    else
        kcatList_merged = kcatList_DLKcat;
    end
     

% STEP 9 Take kcatList and populate edModel.ec.kcat
    ecModel  = selectKcatValue(ecModel, kcatList_merged);

% STEP 10 Get kcat values across isozymes
    ecModel = getKcatAcrossIsozymes(ecModel);

% STEP 11 Calculate molecular weight and store them as model.ec.mw 
    ecModel.ec.mw = zeros(size(ecModel.ec.sequence));
    for j = 1:length(ecModel.ec.mw)
        if strcmp(ecModel.ec.sequence{j},'Sequence unavailable')
            ecModel.ec.mw(j) = 0;
        elseif isnan(ecModel.ec.sequence{j})
            ecModel.ec.mw(j) = 0;
        else
            try
                ecModel.ec.mw(j) = molweight(strrep(ecModel.ec.sequence{j},'X','D'));
            catch 
                fprintf('Protein No. %d does not have molecular weight. Mean molecular weight will be used.\n', j);
                ecModel.ec.mw(j) = 0;
            end
        end
    end
    % Replace the missing protein MWs with avg enzyme weight in this species
    ecModel.ec.mw(ecModel.ec.mw==0) = mean(ecModel.ec.mw(ecModel.ec.mw>0));


% STEP 12 Get standard kcat
    % Fills kcat = 0 (no BRENDA/DLKcat value) with the median kcat of the model,
    % and links reactions without GPR to a 'standard' pseudo-enzyme.
    % The Uniprot database is only used for the standard MW; without it, the
    % median MW of the bin's own proteins is used instead.
    if uniprot_notfound == false
        [ecModel, rxnsMissingGPR, standardMW, standardKcat] = getStandardKcat(ecModel);
    else
        disp("Proteome not found in Uniprot: standard MW taken from the bin's own proteins");
        [ecModel, rxnsMissingGPR, standardMW, standardKcat] = getStandardKcat_noUniprot(ecModel);
    end

% STEP 13 Make spontaneous reactions free of protein cost
    % Keep the 'spontaneous' gene/protein in the ecModel so it can be saved and loaded again,
    % but set kcat = 0 for reactions catalysed only by spontaneous: applyKcatConstraints
    % then gives them zero protein cost, so prot_spontaneous never uses the protein pool.
    % Done after getStandardKcat, otherwise the zero kcats would be refilled with standard kcats.
    spontEnz = strcmpi(ecModel.ec.genes, 'spontaneous');
    if any(spontEnz)
        nEnz = sum(ecModel.ec.rxnEnzMat ~= 0, 2);
        spontEcRxns = any(ecModel.ec.rxnEnzMat(:, spontEnz) ~= 0, 2) & nEnz == 1;
        ecModel.ec.kcat(spontEcRxns) = 0;
        ecModel.ec.source(spontEcRxns) = {'spontaneous'};
        fprintf('%d spontaneous reactions set to zero protein cost.\n', sum(spontEcRxns));
    end

% STEP 14 Apply kcat constraints from ecModel.ec.kcat to ecModel.S
    ecModel = applyKcatConstraints(ecModel);

% STEP 15 Set upper bound of protein pool
    Ptot  = ModelAdapter.params.Ptot;
    f     = ModelAdapter.params.f;
    sigma = ModelAdapter.params.sigma;
    ecModel = setProtPoolSize(ecModel,Ptot,f,sigma);

% STEP 16 Save EC Model
    finalid = strcat(binID, 'Ready_ecGEM.yml');
    model_tmp = fullfile(ecGEMpath, bin, 'models',finalid);
    writeYAMLmodel(ecModel,  model_tmp);

% STEP 17 Run FVA for exchange reactions only 
    % Recommend to use Gurobi solver from RAVEN for FVA 
    % set RAVEN solver as gurobi
    setRavenSolver('gurobi');

    % check the reactions in ecModel ensure exchange reations are completed 
    [ecModel.rxns, constructEquations(ecModel)];
    
    % Extract the exchange reactions
    exchangeRxns = getExchangeRxns(ecModel);
    keep_idx = startsWith(exchangeRxns, "EX_");
    exchangeRxns = exchangeRxns(keep_idx);
    
    % Create matrix to store LB and UB
    n_exc = length(exchangeRxns);
    LB_exc = zeros(n_exc, 1);
    UB_exc = zeros(n_exc, 1); 
    equations= cell(n_exc, 1);

    % FVA
    sensitivity = ;
    sol0 = solveLP(ecModel);
    objFlux = sol0.x(ecModel.c>0);
    objFlux = [objFlux; NaN];
    objFlux = objFlux(1);
    bin = string(bin);

    % Store the growth rate of every bin, including 0 and NA
    tmp = table(sample, bin, objFlux, 'VariableNames', {'sample', 'bin', 'growth'});
    growth_table = [growth_table; tmp];

    if isnan(objFlux) || objFlux <= 0
        disp(strcat("Growth of ", bin, " is 0 or NA, FVA skipped."));
    else
        % run FVA for (sensitivity rate) of the maximun growth rate 
        bioModel = setParam(ecModel,'lb',ecModel.rxns(ecModel.c>0),sensitivity * objFlux);
        bioModel = setParam(bioModel,'ub',ecModel.rxns(ecModel.c>0),objFlux);
        [found, idx] = ismember(exchangeRxns, bioModel.rxns);
        
        for k = 1:length(exchangeRxns)
            iRxn = bioModel.rxns(idx(k));
            iModel = setParam(bioModel,'obj',iRxn,1);
            iSol = solveLP(iModel);
            try
                UB_exc(k) = iSol.x(ismember(iModel.rxns, iRxn));
            catch
                UB_exc(k) = sol0.x(logical(iModel.c));
            end
            
            iModel = setParam(bioModel,'obj',iRxn,-1);
            iSol = solveLP(iModel);
            try
                LB_exc(k) = iSol.x(ismember(iModel.rxns, iRxn));
            catch
                LB_exc(k) = sol0.x(logical(iModel.c));
            end
           
            equations{k} = constructEquations(bioModel, iRxn);
        end

        exchange_results = table(exchangeRxns, equations, LB_exc, UB_exc);
        fva_filename = strcat(sample, ".",binID, 'FVA.txt');
        outputfile = fullfile(ecGEMpath, bin, 'output', fva_filename);
        writetable(exchange_results, outputfile,'Delimiter', '\t', 'QuoteStrings', false, 'FileType', 'text');
    end  

    % Mark the bin as finished: its growth is written last, so a bin that crashed
    % halfway has no growth file and is run again on resume
    fid = fopen(binGrowthFile, 'w');
    fprintf(fid, '%s\t%s\t%.10g\n', sample, bin, objFlux);
    fclose(fid);
end

% Export the growth rate of all bins in this sample
growth_file = fullfile(ecGEMpath, strcat(sample, "_growth.txt"));
writetable(growth_table, growth_file,'Delimiter', '\t', 'QuoteStrings', false, 'FileType', 'text');
