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

minGlucose = 0.1; %Need some glucose to grow in some models.

if ~exist('linPath', 'vars')
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

O2Upt = rxnIdxFromNames(mRev, OxygenRxns); 
mRev.lb(O2Upt) = -2; % low default found in iMM904

pyrUpt = rxnIdxFromNames(mRev, PyruvateRxns); 
mRev.lb(pyrUpt) = -1; % set for convertToIrreversible

glucUpt = rxnIdxFromNames(mRev, GlucoseRxns); 
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
pdc = zeros(length(pyrRxns), nIter);
eth = zeros(length(ethRxns) + 1, nIter);
qMat = zeros(nIter, pathLen - 1);
fluxMat = zeros(nRxns, nIter);
biomvec = zeros(1, nIter);

pathInfo = [];
for i =1:nIter
  m.ub(itermet) = 1 * (initF + (i-1)*step);
  disp(m.ub(itermet));
  solt = optimizeCbModel(m, 'max', 'one');
  %pdcratio = solt.x(pyrRxns)/sum(solt.x(pyrRxns));
  %pdc(:, i) = pdcratio';
  [qVec, pathInfo] = linPathFlux(mIrr, linPath, solt.x);
  pdc(:, i) = solt.x(pyrRxns); %easy to calculate ratio later.
  ethsum = sum(solt.x(ethRxns));
  eth(:, i) = [solt.x(ethRxns)' ethsum]';
  biomvec(i) = solt.f;
end

save(['PyruvateAll_' m.rxnNames{itermet} '_' num2str(nIter) '.mat'], ...
      'qMat', 'fluxMat', 'pathInfo');


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

