function [mutMatBlock, ActualBlkSizes] = sortRandSampleMutants(model, grRateKOInc, WTensInc, grRateKODec, WTensDec, mutMatMax, Nint, minfit, maxfit, BlockSize)
%Pick random values within intervals
%Nint_bnd_i = minfit+i*(maxfit-minfit)/Nint
mutMatBlock = zeros(Nint*BlockSize,3);
ActualBlkSizes = zeros(1,Nint);
incSize = (maxfit-minfit)/Nint;

allRxns = 1:length(model.rxns);
geneRxns = allRxns(any(model.rxnGeneMat,2));


grTens = [grRateKOInc grRateKODec];
WTens = [WTensInc WTensDec];
grsz = size(grTens)


%There are two types of dubious mutants; those with weights in excess of the optimal
%mutant weight, and those with the opposite sign. Take care of (at least the first)
%of these here by editing grTens.
for i = 1:length(mutMatMax)
  for j = 1:grsz(2)
    if (grTens(mutMatMax(i,3),j) > minfit) && ( (mutMatMax(i,2) > 0 && (WTens(mutMatMax(i,3),j) > mutMatMax(i,2)) ) ...
       || (mutMatMax(i,2) < 0 && (WTens(mutMatMax(i,3),j) < mutMatMax(i,2)) ) )
      grTens(mutMatMax(i,3),j) = 0;      
    end
    %Also take care of opposite sign case:
    if WTens(mutMatMax(i,3),j) * mutMatMax(i,2) < 0
      grTens(mutMatMax(i,3),j) = 0;
    end
  end
end


%Need to make a big tensor with all the data in it, then sort by fitness.
%Sort the fitnesses with their idx, then create the tensor based on this sorting,
%since matlab can't sort a tensor in this way.
BigTens = zeros(grsz(1),grsz(2),3);
fitIdx = zeros(grsz(2),2);
for i = 1:grsz(1)
  fitIdx(:,1) = grTens(i,:);
  fitIdx(:,2) = 1:grsz(2);
  fitIdx = sortrows(fitIdx,1);
  for j = 1:grsz(2)
    BigTens(i,j,1) = fitIdx(j,1);
    BigTens(i,j,2) = WTens(i,fitIdx(j,2)); 
    BigTens(i,j,3) = i; 
  end
end

%Next, for each interval and each gene, find the minimum and maximum fitness indices.
IntMin = zeros(grsz(1),Nint);
IntMax = zeros(grsz(1),Nint);
for i = 1:grsz(1)
  lastMin = 1;
  lastMax = 1;
  for j = 1:Nint
    bndNintLow = minfit + (j-1)*incSize;
    bndNintHigh = minfit + j*incSize;
    for k=lastMin:grsz(2)
      if BigTens(i,k,1) >= bndNintLow
        IntMin(i,j) = k;
        lastMin = k;
        break
      end
    end 
    for k=lastMax:grsz(2)
      if BigTens(i,k,1) >= bndNintHigh
        if (BigTens(i,k-1,1) < bndNintHigh) && (BigTens(i,k-1,1) > (bndNintHigh-incSize)) 
          IntMax(i,j) = k-1;
          lastMax = k-1;
        end
        break;
      elseif k==grsz(2) && BigTens(i,k,1) > bndNintLow
        IntMax(i,j) = k;
      end
    end
  end
end

%Now, randomly sample (without resampling) genes that match the interval criteria.
%At the same time, randomly pick within the gene's interval a mutant, and create mutMatBlock.
%rvec = randNoResample(rmax, rlen, excludes)
for i = 1:Nint
  %get suitable genes; it makes sense I'm pretty sure:
  candidateGenes = find(((IntMax(:,i) - IntMin(:,i)+1) .* IntMin(:,i) > 0));
  candidateGenes = intersect(candidateGenes, geneRxns); %Should optimize this above
  size(candidateGenes)
  rvec = randNoResample(length(candidateGenes), BlockSize);
  sampledGenes = candidateGenes(rvec);
  ActualBlkSizes(i) = length(sampledGenes); 
  for j = 1:length(sampledGenes)
    rrange = IntMax(sampledGenes(j),i)-IntMin(sampledGenes(j),i)+1;
    randmut = IntMin(sampledGenes(j),i) + randi(rrange) - 1;
    mutMatBlock((i-1)*BlockSize+j,:) = BigTens(sampledGenes(j),randmut,:);
  end
end


%Since in each block the fitnesses may not be sorted, sort them here.
mutMatBlock = sortrows(mutMatBlock,1);
%Remove 0-value entries
notzvals = find(mutMatBlock(:,1)~=0);
mutMatBlock = mutMatBlock(notzvals,:);
%Consider how to take an average of these to make a nicer heatmap.

