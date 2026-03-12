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

%% Set model adapter
adapter = ModelAdapterManager.setDefault('/Users/zhangtianyue/Desktop/genome/GSM/M/GECKO/M_newModelGEM/M_newModelGEMAdapter.m');

%% Load conventional GEM model
model = loadConventionalGEM();

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

%% Load EC model at 29°C condition
ecModel = loadEcModel('Mtaihu98_light_ecModel.yml');

%% Import new kcat values
% kcat data file contains a single column without header
newKcat = readmatrix('M15_kcat.tsv', 'FileType', 'text');

%% Replace existing kcat values in EC model
ecModel.ec.kcat = newKcat;

%% Save model for manual adjustment of exchange reaction bounds
saveEcModel(ecModel,'Mtaihu98_T15_ecModel.yml' ); % Update filename if necessary

%% Reload EC model after modifications
ecModel = loadEcModel('Mtaihu98_T15_ecModel.yml');

%% Apply kcat constraints to stoichiometric matrix
ecModel = applyKcatConstraints(ecModel);

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
writetable(fluxTable, 'M15_GECKO_light_Flux_samePro.csv');
