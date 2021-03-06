function [epiMats smfVecs] = multipleEpistasisMatsPart(model, modelNew, atime, WTflux, grRateKOInc, WTensInc, grRateKODec, WTensDec, mutMatMax, Nint, minfit, maxfit, BlockSize, epiIter)
partWidth = 100;

[mutMatBlock, ActualBlkSizes] = sortRandSampleMutants(model, grRateKOInc, WTensInc, grRateKODec, WTensDec, mutMatMax, Nint, minfit, maxfit, BlockSize);
size(mutMatBlock)
mutMatBlock = squeeze(mutMatBlock(1:partWidth,:));
[grMat, epiMult, smfAdj] = doubleRxnMutationObjWithin(model, modelNew, atime, mutMatBlock, WTflux, ones(1,length(model.lb)), true);

ngmuts = length(smfAdj);

epiMats = zeros(epiIter,ngmuts,ngmuts);
smfVecs = zeros(epiIter,ngmuts);
for i = 1:(epiIter-1)
  disp(strcat('Finished ',num2str(i)));
  epiMats(i,:,:) = epiMult;
  smfVecs(i,:) = smfAdj;
  save('epiMatsPart.mat','epiMats','smfVecs');
  [mutMatBlock, ActualBlkSizes] = sortRandSampleMutants(model, grRateKOInc, WTensInc, grRateKODec, WTensDec, mutMatMax, Nint, minfit, maxfit, BlockSize);
  [grMat, epiMult, smfAdj] = doubleRxnMutationObjWithin(model, modelNew, atime, squeeze(mutMatBlock(1:partWidth,:)), WTflux, ones(1,length(model.lb)), true);
  

end

epiMats(i+1,:,:) = epiMult;
smfVecs(i+1,:) = smfAdj;
save('epiMatsPart.mat','epiMats','smfVecs');