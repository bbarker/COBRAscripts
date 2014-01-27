function [pdc, biomvec] = SSPyruvateAll(model, itermet, initF, finalF, step)
%Assumes an irreversible model!
%O2 uptake: rxn 882
%glucose  : rxn 794
%pyruvate : rxn 918
%PYRDC    : rxn 1958
pyrRxns = [140 144 238 1762 1770 1958]; %Exclude 1962, a transport rxn.
mm = model;
nIter = (finalF-initF)/step + 1;
mm.lb(itermet) = initF;
pdc = zeros(length(pyrRxns),nIter);
biomvec = zeros(1,nIter);

for i =1:nIter
  mm.ub(itermet) = 1 * (initF + (i-1)*step);
  disp(mm.ub(itermet));
  solt = optimizeCbModel(mm,'max',true,true);
  pdcratio = solt.x(pyrRxns)/sum(solt.x(pyrRxns));
  pdc(:,i) = pdcratio';
  biomvec(i) = solt.f;
end





