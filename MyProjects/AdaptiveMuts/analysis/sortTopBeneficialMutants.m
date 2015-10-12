function mutMat = sortTopBeneficialMutants(modelOld,model,method,WTflux,objMultInit,objMult,dlvl,FVAmin,FVAmax,grRateKOTens,SPvec)

%Only interested in enzyme-associated reactions.                                                                              
if nargin < 11
  SPvec=ShadowPricesHeuristic(modelOld,WTflux);
end

grsz = size(grRateKOTens);
%nmut = grsz(2)-1;
 
allRxns = 1:length(model.rxns);
geneRxns = allRxns(any(model.rxnGeneMat,2));

[solAda,solutionWT,totalFluxDiff,solStatus] = MOMATest(modelOld,model,'max',false,true,objMultInit,WTflux,SPvec);

%Remove WT values:
WTidx = grsz(2)/2;
size(grRateKOTens(:,1:(WTidx-1)))
size(grRateKOTens(:,(WTidx+1):grsz(2)))
size(grRateKOTens)

maxMVal = max(grRateKOTens');
maxBounds = zeros(1,grsz(1));
maxGeneNum = zeros(1,grsz(1));
mutMat = zeros(3,grsz(1));

%Actually don't return 3 arguments, use sort() instead
parfor i = 1:grsz(1)
  Idx = find(squeeze(grRateKOTens(i,:))==maxMVal(i));
  mval1 = Idx(1);
  if mval1 == ceil((grsz(2)-1)/2)
    mv1part = -1;
  elseif mval1 < (grsz(2)-1)/2
    mv1part = (mval1-1)/floor((grsz(2)-1)/2);
  elseif mval1 > (grsz(2)-1)/2      
    mv1part = (mval1-ceil((grsz(2)-1)/2))/floor((grsz(2)-1)/2);  
  end

  modelTmp = model;
  %if mv1 > 0
  if mval1 > 0
      if mv1part < 0                    
	modelTmp.ub(i) = solAda.x(i);
	modelTmp.lb(i) = solAda.x(i);
      elseif mval1 < (grsz(2)-1)/2
	modelTmp.ub(i) = solAda.x(i)+(1-mv1part)*(FVAmin(i)-solAda.x(i));
	modelTmp.lb(i) = solAda.x(i)+(1-mv1part)*(FVAmin(i)-solAda.x(i));                    
      else
	modelTmp.ub(i) = solAda.x(i)+(mv1part)*(FVAmax(i)-solAda.x(i));
	modelTmp.lb(i) = solAda.x(i)+(mv1part)*(FVAmax(i)-solAda.x(i));
      end
  end
  %if mv1 == 0
  if mval1 > (grsz(2)-1) %just stick deletions on the end
    modelTmp.ub(i) = solAda.x(i)*dlvl;
    modelTmp.lb(i) = solAda.x(i)*dlvl;
  end

  [solKO,solutionWT,totalFluxDiff,solStatus] = MOMATest(modelOld,modelTmp,'max',false,true,objMult,WTflux,SPvec);
  if abs((solKO.f - grRateKOTens(i,mval1))/solAda.f)> 1e-3
    error('Gene %d and mut %d did not match prior data', i, mval1);
  end

  maxMVal(i) = solKO.f;
  maxBounds(i) = modelTmp.lb(i);  
  maxGeneNum(i)=i; 
end


mutMat = [maxMVal; maxBounds; maxGeneNum];
mutMat = mutMat(:,geneRxns);
mutMat = sortrows(mutMat',1);
%grRateKOTens = [grRateKOTens(:,1:(WTidx-1)) grRateKOTens(:,(WTidx+1):grsz(2))];

