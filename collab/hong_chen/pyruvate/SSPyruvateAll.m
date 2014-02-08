function [qMat, fluxMat] = SSPyruvateAll(m, itermet, initF, finalF, step, linPath)
% Assumes an irreversible model!

%
% This script records the fraction of pyruvate converted to acetaldehyde
% (PDC) and then in to ethanol in an FBA simulation where increasing
% levels of glucose (or pyruvate, or some other metabolite) are allowed
% into the system.
%

%
% Glucose Uptake Rxns
% iMM904 : 'D Glucose exchange backward'
% Yeast7 : 'D-glucose exchange backward'
GlucoseRxns = {'D Glucose exchange backward', 'D-glucose exchange backward'};

%
% Pyruvate Uptake Rxn
% iMM904 :
% Yeast7 :
%

%
%Need to set major uptake levels to be the same for various models
%
% Carbon source (Glucose, Pyruvate), O2;  phosphorous, ammonia ?
%
%
m.ub(...) = ...

% itermet : usually glucose (rxn = 727), pyruvate (rxn = 797)

nRxns = length(m.rxns);
if itermet == 727
    m.ub(797) = 0;
end
if itermet == 797
    m.ub(727) = 0.1; %need some glucose
end
if ~exist('linPath', 'vars')
    linPath = {'Pyruvate', 'Acetaldehyde', 'Ethanol'}; %iMM904
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





