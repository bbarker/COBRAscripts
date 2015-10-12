function [epiMats smfVecs] = multipleEpistasisMats(model, modelNew, atime, WTflux, grRateKOInc, WTensInc, grRateKODec, WTensDec, mutMatMax, Nint, minfit, maxfit, BlockSize, epiIter)

[mutMatBlock, ActualBlkSizes] = sortRandSampleMutants(model, grRateKOInc, WTensInc, grRateKODec, WTensDec, mutMatMax, Nint, minfit, maxfit, BlockSize);
[grMat, epiMult, smfAdj] = doubleRxnMutationObjWithin(model, modelNew, atime, mutMatBlock, WTflux, ones(1,length(model.lb)), true);

ngmuts = length(smfAdj);

epiMats = zeros(epiIter,ngmuts,ngmuts);
smfVecs = zeros(epiIter,ngmuts);
for i = 1:(epiIter-1)
  disp(strcat('Finished ',num2str(i)));
  epiMats(i,:,:) = epiMult;
  smfVecs(i,:) = smfAdj;
  save('epiMats.mat','epiMats','smfVecs');
  [mutMatBlock, ActualBlkSizes] = sortRandSampleMutants(model, grRateKOInc, WTensInc, grRateKODec, WTensDec, mutMatMax, Nint, minfit, maxfit, BlockSize);
  [grMat, epiMult, smfAdj] = doubleRxnMutationObjWithin(model, modelNew, atime, mutMatBlock, WTflux, ones(1,length(model.lb)), true);
  

end

epiMats(i+1,:,:) = epiMult;
smfVecs(i+1,:) = smfAdj;
save('epiMats.mat','epiMats','smfVecs');