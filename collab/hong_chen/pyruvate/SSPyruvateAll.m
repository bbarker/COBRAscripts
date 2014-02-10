function [qMat, fluxMat] = SSPyruvateAll(mRev, itermet, initF, finalF, step, linPath)
% Assumes a reversible model! (which is then later converted to irrev)

%
% This script records the fraction of pyruvate converted to acetaldehyde
% (PDC) and then in to ethanol in an FBA simulation where increasing
% levels of glucose (or pyruvate, or some other metabolite) are allowed
% into the system.
%
% itermet : usually a carbon source (e.g. 'glucose' or 'pyruvate')
%

minGlucose = 0.5; %Need some glucose to grow in some models.

if ~exist('linPath', 'var')
    linPath = {'Pyruvate', 'Acetaldehyde', 'Ethanol'}; %iMM904
end

% Glucose Uptake Rxns
% iMM904 : 'D Glucose exchange'
% Yeast7 : 'D-glucose exchange'
GlucoseRxns = {'D Glucose exchange', 'D-glucose exchange'};
GlucoseRxnsBack = rNameBack(GlucoseRxns);
%
% Pyruvate Uptake Rxn
% iMM904 : 'Pyruvate exchange'
% Yeast7 : 'pyruvate exchange'
%
PyruvateRxns = {'Pyruvate exchange', 'pyruvate exchange'};
PyruvateRxnsBack = rNameBack(PyruvateRxns);

%
% Oxygen Uptake Rxn
% iMM904 : 'O2 exchange'
% Yeast7 : 'oxygen exchange'
%
OxygenRxns = {'O2 exchange', 'oxygen exchange'};
OxygenRxnsBack = rNameBack(OxygenRxns);

%
%Need to set major uptake levels to be the same for various models
%
% Carbon source (Pyruvate uptake), O2;  phosphorous, ammonia ?
%
%

O2UptRev = rxnIdxFromNames(mRev, OxygenRxns); 
mRev.lb(O2UptRev) = -2; % low default found in iMM904
%
pyrUptRev = rxnIdxFromNames(mRev, PyruvateRxns); 
mRev.lb(pyrUptRev) = -1; % set for convertToIrreversible
%
glucUptRev = rxnIdxFromNames(mRev, GlucoseRxns); 
%no need to set it here

%
%
% Now in irreversible land
%
%
m = convertToIrreversible(mRev);
nRxns = length(m.rxns);

O2Upt = rxnIdxFromNames(m, OxygenRxnsBack); 
pyrUpt = rxnIdxFromNames(m, PyruvateRxnsBack); 
glucUpt = rxnIdxFromNames(m, GlucoseRxnsBack); 

%Not enough options to worry about using a map for now.
if strcmp(itermet, 'pyruvate')
    itermet = pyrUpt;
    m.ub(glucUpt) = minGlucose; %need some glucose
elseif strcmp(itermet, 'glucose')
    itermet = glucUpt;
    m.ub(pyrUpt) = 0;
end

nIter = floor((finalF-initF)/step) + 1;
pathLen = length(linPath);
m.lb(itermet) = initF;
qMat = zeros(nIter, pathLen - 1);
fluxMat = zeros(nRxns, nIter);

pathInfo = [];
parfor i = 1:nIter
    mt = m;
    mt.ub(itermet) = 1 * (initF + (i-1)*step);
    disp(mt.ub(itermet));
    solt = optimizeCbModel(mt, 'max', 'one');
    [qVec, pathInfo] = linPathFlux(mt, linPath, solt.x);
    qMat(i, :) = qVec;
    fluxMat(:, i) = solt.x;
end

modName = strrep(strrep(strrep(strrep(num2str(m.description), ' ', ''), ...
          '.', ''), 'xml', ''), '_', '');

save(['PyruvateAll_' modName '_' m.rxnNames{itermet} '_' ...
      num2str(nIter) '.mat'], 'qMat', 'fluxMat', 'pathInfo');


end % of SSPyruvateAll

function rxn = rxnIdxFromNames(m, rNames)
% Given a list of rxnNames for the same reaction in various
% models, return the reaction index for the first found rNames
% in m.rxnNames.
[~, rIdx] = ismember(rNames, m.rxnNames);
rIdxNotZero = find(rIdx);
rxn = rIdx(rIdxNotZero(1));
end % of rxnIdxFromNames

function rNamesIrr = rNameBack(rNames)
rNamesIrr = cellfun(@(x) strcat(x, ' backward'), rNames, ... 
    'UniformOutput', false);
end % of rNameBack

