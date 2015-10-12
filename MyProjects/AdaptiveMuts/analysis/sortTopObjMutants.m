function mutMat = sortTopObjMutants(model, grRateKOInc, VTensInc, WTensInc, grRateKODec, VTensDec, WTensDec)

allRxns = 1:length(model.rxns);
geneRxns = allRxns(any(model.rxnGeneMat,2));


grTens = [grRateKOInc grRateKODec];
%VTens = [VTensInc VTensDec];
WTens = [WTensInc WTensDec];
grsz = size(grTens);

%mutMat = zeros(3,2*grsz(1)*grsz(2));


disp(size([VTensInc VTensDec]));
maxMVal = max(grTens');
size(maxMVal)
maxWeight = zeros(1,grsz(1));
maxGeneNum = zeros(1,grsz(1));
mutMat = zeros(3,grsz(1));

parfor i = 1:grsz(1)
  Idx1 = find(squeeze(grTens(i,:))==maxMVal(i));
  %WvalsMap = [1./(WTensInc-1) WTensDec];
  WvalsMap = [WTensInc abs(WTensDec)]
  minmaxval = max(WvalsMap(Idx1));
  %Intuitively, it makes sense to choose the lowest weight (at least to me)
  %but data suggests this may not be the best option for combinations;
  %try the max
  %minmaxval = min(WvalsMap(Idx1));
  Idx2 = find(WvalsMap==minmaxval);
  Idx = intersect(Idx1, Idx2);
  % maxMVal(i) = solKO.f;
  maxWeight(i) = WTens(i,Idx(1));
  maxGeneNum(i)=i;   
end


%WTensInc = flat([WTensInc; WTensDec]);
%mutMat = zeros(3,2*grsz(1)*grsz(2));
mutMat = [maxMVal; maxWeight; maxGeneNum];
mutMat = mutMat(:,geneRxns);
mutMat = sortrows(mutMat',1);
%grRateKOTens = [grRateKOTens(:,1:(WTidx-1)) grRateKOTens(:,(WTidx+1):grsz(2))];

