function [mutDat, numFound] = randomSingleBeneMuts(model, InitFlux, minFlux, maxFlux, sthresh, initR, sampType, fracInBnd)
rmpath('/home/brandon/FBA/mosek/6/toolbox/r2009b/');
doDel = true;
objWeight=0; %!!!!
maxIter = 50;
%!!Remember to turn on maxIterTmp += ...
zthresh = 1e-8;
nDesired = 50;
growthRxn = find(model.c);

if nargin<6
  initR = 1;
end

if nargin<7
  sampType = 'uniform'
end

if nargin<8
  fracInBnd = 0.95
end

%Use == for mutation restrictions.
%
%ToFlux: optimal not needed now, may use later in biased sampling
%minFlux, maxFlux: output from (loopless) FVA
%maxIter: how many epistasis values to attempt.
%mutDat has the following columns (fx is relative):
%  (fitness, flux restriction, reaction #)



nRxns = length(InitFlux);
geneRxns = any(model.rxnGeneMat,2);
geneRxns = find(geneRxns);
nGRxns = length(geneRxns);
numFound = zeros(1, nGRxns);

mutDat = zeros(nGRxns,nDesired,3);
mutDatDel = zeros(nGRxns,nDesired,3);
tmpvec = zeros(nDesired,3);
rvec = zeros(nDesired,1);
zvec = zeros(1,nRxns);

for i = initR:nGRxns
  iter = 0;
  filled = 0;
  del_filled = 0;
  rxnNum1 = geneRxns(i);
  FVAwidth = maxFlux(rxnNum1) - minFlux(rxnNum1);

  disp(['Starting Rxn #' num2str(i)]); 
  maxIterTmp = maxIter;
  if FVAwidth > zthresh
    sigma = 0;
    if strcmp(sampType,'truncnorm')
      sigma = findSigma(InitFlux(rxnNum1), minFlux(rxnNum1), maxFlux(rxnNum1), fracInBnd);
    end
    while (filled < nDesired) & (iter < maxIterTmp)
      tmpvec = zeros(nDesired,3);
      %Sampling type: edit the following for desired distribution.
      switch sampType
        case 'uniform'
      %Uniform:
          rvec = minFlux(rxnNum1) + FVAwidth*rand(nDesired,1);
	case 'truncnorm'
      %Truncated normal. First find sigma, then compute rvec:
          rvec = simpleTruncatedNorm(sigma, minFlux(rxnNum1), maxFlux(rxnNum1), nDesired, InitFlux(rxnNum1));
	otherwise 
	  disp('error: no valid sample type');
      end
      parfor j = 1:nDesired
        modelMut = model;
	modelMut.lb(rxnNum1) = rvec(j);
	modelMut.ub(rxnNum1) = rvec(j);
        [solutionDel,solutionWT] = linearMOMAOneShot(modelMut,InitFlux,objWeight,0,1,zvec);
	tmpvec(j,:) = [solutionDel.f/InitFlux(growthRxn) rvec(j) rxnNum1];       
      end
      iter = iter + 1;
      for j = 1:nDesired
	if tmpvec(j,1) > 1+sthresh
	  filled = filled+1;
	  maxIterTmp = maxIterTmp + 5*maxIter;
	  if filled<=nDesired
	    mutDat(i,filled,:) = tmpvec(j,:);
	  end
	elseif (tmpvec(j,1) <= 1-sthresh) & tmpvec(j,1) > zthresh
	  del_filled = del_filled+1;
	  if del_filled<=nDesired
	    mutDatDel(i,del_filled,:) = tmpvec(j,:);
          end
	end
      end
    end
  end
  numFound(i) = min(filled,nDesired);
  numDelFound(i) = min(del_filled,nDesired);
  lasti = i;
  save('randomSingleBeneMuts_tmp.mat', 'mutDat', 'mutDatDel', 'numFound','numDelFound','lasti');
end


%Perform Sorting By Fitness!!!!
if doDel
  mutDat = [mutDat mutDatDel];
  numFound = numFound + numDelFound;
end
parfor i = 1:nGRxns
  mutDat(i,:) = sortrows(squeeze(mutDat(i,:)),1);
end
