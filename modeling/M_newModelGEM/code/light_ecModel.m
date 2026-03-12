%% Clear workspace
clear all

%% Add required paths
addpath(genpath('/Users/zhangtianyue/Desktop/GECKO/Tool/GECKO-3.2.3'));
addpath(genpath('/Users/zhangtianyue/Desktop/GECKO/Tool/RAVEN-2.10.3'));  % Add RAVEN toolbox if used
savepath;  % Save path configuration for future MATLAB sessions

%% Set solver
changeCobraSolver('gurobi', 'all');

%% Install GECKO toolbox
cd('/Users/zhangtianyue/Desktop/GECKO/Tool/GECKO-3.2.3')
GECKOInstaller.install

%% Initialize GECKO project (required for new models)
% After running, update adapter path and configuration if necessary
startGECKOproject();

%% Set model adapter
adapter = ModelAdapterManager.setDefault('/Users/zhangtianyue/Desktop/genome/GSM/M/GECKO/M_newModelGEM/M_newModelGEMAdapter.m');

%% Load conventional GEM model
model = loadConventionalGEM();

%% Preserve original metabolite compartment information
old_metComps = model.metComps;

%% Convert model to RAVEN-compatible structure
model = ravenCobraWrapper(model);

%% Manually import GPR information
model.grRules = readlines('grRules.tsv');
model.genes = readlines('genes.tsv');
model.eccodes = readlines('eccodes.tsv');

%% Convert data types to cell arrays
model.grRules = cellstr(model.grRules);
model.genes = cellstr(model.genes);
model.eccodes = cellstr(model.eccodes);

%% Rebuild reaction-gene mapping matrix
model = buildRxnGeneMat(model);

%% Change directory to EC model project (optional)
cd('/Users/zhangtianyue/Desktop/genome/GSM/M/GECKO/ecM_Tahihu98GEM')
addpath('/Users/zhangtianyue/Desktop/genome/GSM/M/GECKO/ecM_Tahihu98GEM')
savepath

%% Restore original metabolite compartment information
model.metComps = old_metComps;

%% Create EC model from GEM
[ecModel, noUniprot] = makeEcModel(model, true);

%% Assign model identifier
ecModel.id = 'Mtaihu98';

%% Import metabolite SMILES information
ecModel.metSmiles = readlines('metSmiles.tsv');
ecModel.metSmiles = cellstr(ecModel.metSmiles);

%% Import metabolite names
ecModel.metNames = readlines('metNames.tsv');
ecModel.metNames = cellstr(ecModel.metNames);

%% Retrieve EC data from BRENDA
ecModel = getECfromGEM(ecModel);
kcatList_fuzzy = fuzzyKcatMatching(ecModel);

%% Filter unrealistically high kcat values
% Values greater than 1000 are set to zero for DLKcat prediction
kcatList_fuzzy.kcats(kcatList_fuzzy.kcats > 1000) = 0;

%% Predict kcat values using DLKcat
ecModel = findMetSmiles(ecModel);

%% Generate DLKcat input files
writeDLKcatInput(ecModel,[],[],[],[],true);

%% Run DLKcat prediction
% Requires local DLKcat installation and configuration
runDLKcat();

%% Load DLKcat output results
kcatList_DLKcat = readDLKcatOutput(ecModel);

%% Merge BRENDA-derived and DLKcat-derived kcat values
kcatList_merged = mergeDLKcatAndFuzzyKcats(kcatList_DLKcat, kcatList_fuzzy);

%% Implement selected kcat values into EC model
ecModel = selectKcatValue(ecModel, kcatList_merged);

%% Retrieve standard kcat and molecular weight values (optional)
[ecModel, rxnsMissingGPR, standardMW, standardKcat] = getStandardKcat(ecModel);

%% Apply kcat constraints to stoichiometric matrix
ecModel = applyKcatConstraints(ecModel);

%% Save EC model
saveEcModel(ecModel, 'Mtaihu98_light_ecModel.yml');

%% Load EC model for further analysis
ecModel = loadEcModel('Mtaihu98_light_ecModel.yml');

%% Set protein pool upper bound
% Default upper limit: 0.5 g/gDW (model-defined)
ecModel = setProtPoolSize(ecModel, 0.5, 0.5, 0.5);

%% Solve linear programming problem
sol = solveLP(ecModel);
printFluxes(ecModel, sol.x);

%% Define objective function and optimize growth
ecModel = setParam(ecModel,'obj',adapter.params.bioRxn,1);
sol = solveLP(ecModel);
fprintf('Growth rate reached: %g /hour.\n', sol.f);

%% Constrain growth to 99% of optimal value
ecModel = setParam(ecModel,'lb',adapter.params.bioRxn,0.99*sol.f);

%% Minimize protein pool usage
ecModel = setParam(ecModel,'obj','prot_pool_exchange',1);
sol = solveLP(ecModel);
fprintf('Minimum protein pool usage: %g mg/gDCW.\n', sol.f);

%% Map fluxes back to conventional GEM model
[mappedFlux, enzUsageFlux, usageEnz] = mapRxnsToConv(ecModel, model, sol.x);

%% Export reaction IDs and flux values to table
fluxTable = table(model.rxns, mappedFlux, 'VariableNames', {'ReactionID', 'Flux'});

%% Save results as CSV file
writetable(fluxTable, 'M29_GECKO_light_Flux.csv');
