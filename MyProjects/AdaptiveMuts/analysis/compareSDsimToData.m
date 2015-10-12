function [minFitDiff FitDiffsSum FitDiffsCnt RelFitDiffs] = compareSDsimToData(model, fitVec, modelGeneIdx, RelativeFitness)
%FitDiffsSum minimize the summed fitness difference from experimental data.
%FitDiffsCnt minimize the number of mutant fitness differences < perDiff (fraction).
%
%call makeFitVec(model,fitCell,condCol) first.
perDiff = 0.1

numGShared = sum(modelGeneIdx)
sharedIdx = find(modelGeneIdx);
FitDim = size(RelativeFitness);
RelFitDiffs = zeros(FitDim(1), numGShared); 
minFitDiff = zeros(1,FitDim(1));
%RelativeFitness = RelativeFitness(:,sharedIdx);
FitDiffsSum = zeros(1,FitDim(1));
FitDiffsCnt = zeros(1,FitDim(1));
for i = 1:FitDim(1)
  for j = 1:numGShared
    RelFitDiffs(i,j) = abs(RelativeFitness(i,sharedIdx(j))/fitVec(sharedIdx(j))-1); 
  end %Why are NaNs here?
  baddies = isinf(RelFitDiffs(i,:)) | isnan(RelFitDiffs(i,:));
  FitDiffsSum(i) = sum(RelFitDiffs(i,~baddies));
  FitDiffsCnt(i) = sum(RelFitDiffs(i,~baddies)<perDiff);
end

for i = 1:numGShared
  baddies = isinf(RelFitDiffs(:,i)) | isnan(RelFitDiffs(:,i));
  minFitDiff(i) = min(RelFitDiffs(~baddies,i));
end

