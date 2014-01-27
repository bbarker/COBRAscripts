function [pdc, biomvec, eth] = SSPyruvateAll(model, itermet, initF, finalF, step)
%Assumes an irreversible model!
%O2 uptake: rxn 
%glucose  : rxn 
%pyruvate : rxn 
%PYRDC    : rxn 
%ethanol  : rxns   [143 261 262]
%
% itermet : usually glucose (rxn = 727), pyruvate (rxn = 797)
mm = model;
pyrRxns = [238 1609 1805 1806 140 144 240 603 1617]; %Exclude transport rxns.
                                   %^ These are for iMM904.
ethRxns = [143 261 262];
mm.ub(264) = 0; % 264 ends up being a duplicate of 261 
if itermet == 727
    mm.ub(797) = 0;
end
if itermet == 797
    mm.ub(727) = 0.1; %need some glucose
end
nIter = floor((finalF-initF)/step) + 1;
mm.lb(itermet) = initF;
pdc = zeros(length(pyrRxns), nIter);
eth = zeros(length(ethRxns) + 1, nIter);
biomvec = zeros(1, nIter);

for i =1:nIter
  mm.ub(itermet) = 1 * (initF + (i-1)*step);
  disp(mm.ub(itermet));
  solt = optimizeCbModel(mm, 'max', 'one');
  pdcratio = solt.x(pyrRxns)/sum(solt.x(pyrRxns));
  pdc(:, i) = pdcratio';
  ethsum = sum(solt.x(ethRxns));
  eth(:, i) = [solt.x(ethRxns)' ethsum]';
  biomvec(i) = solt.f;
end

save(['PyruvateAll_' mm.rxnNames{itermet} '_' num2str(nIter) '.mat'], ...
    'pyrRxns', 'pdc', 'biomvec', 'ethRxns', 'eth');





