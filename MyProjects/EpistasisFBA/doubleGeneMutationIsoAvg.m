function grRateKOTens = doubleGeneMutationIsoAvg(model,method, fitvals1, fitvals2, WTflux, dlvl, geneList1,geneList2,enzConstraint)
%doubleGeneMutation Performs double gene deletion analysis using FBA, MOMA,
%or linear MOMA
%
  % grRateKOTens =  doubleGeneDeletion(model,method,geneList1,geneList2,verbFlag)
%
% model     Model structure including gene-reaction associations
% method    Either 'FBA' (default) 'MOMA', or 'lMOMA' (optional)
% geneList1 List of query genes to be deleted (optional, default is all genes)
% geneList2 List of target genes to be deleted (optional, default is
% geneList1)
% verbFlag  Verbose output (optional,default false)
%
% type
% grRateKOTens  Double deletion strain growth rates (1/h)
% grRateWT  Wild type growth rate (1/h)
%
% Markus Herrgard 8/8/06

%Assumptions: 1) All genes listed are in the model
% 2) the objective is the biomass or growth rate, but here we always call it growth rate

if (nargin < 2)
    method = 'FBA';
end
if (nargin < 7)
    geneList1 = model.genes;
    differentSetsFlag = false;
else
    if (isempty(geneList1))
        geneList1 = model.genes;
    end
end

if (nargin < 8)
    geneList2 = geneList1;
    differentSetsFlag = false;
else
    if (isempty(geneList2))
        geneList2 = geneList1;
        differentSetsFlag = false;
    else
        differentSetsFlag = true;
    end
end
if (nargin < 9)
    enzConstraint=1;
end

%fitvals=[0 0.25 0.5 0.75 1];
%fitvals=[0 0.5 1 1.5 2];
%fitvals1=[0 1];
%fitvals2=[0 1];

nGenes = length(model.genes);
%solTruMax = optimizeCbModel(model,'max',true,true);
%grRateMax = solTruMax.f;
if enzConstraint ~= 1
  disp('blah')
  model = constrainEnzFluxGlobal(model, enzConstraint);
end

solMax = optimizeCbModel(model,'max');
grRateWT = solMax.f;
if nargin < 5
  solMax = optimizeCbModel(model,'max',true,true);
  WTflux = solMax.f;
%else
%  disp('Using precomputed flux vector');
end
grRateWT = solMax.f;



nDelGenes1 = length(geneList1);
nDelGenes2 = length(geneList2);

%grRateKOTens = ones(nDelGenes1,nDelGenes2, 5*5)*grRateWT;
grRateKOTens = zeros(nDelGenes1, nDelGenes2, length(fitvals1), length(fitvals2));

if (differentSetsFlag)
    nTotalPairs = nDelGenes1*nDelGenes2;
else
    nTotalPairs = nDelGenes1*(nDelGenes1-1)/2;
end

% Run double deletion analysis


for mval1 = 1:length(fitvals1)
   mv1 = fitvals1(mval1);
   for mval2 = 1:length(fitvals2) 
      mv2 = fitvals2(mval2);
      if mv2 == 1 && mv1 == 1
         continue
      end
      grRateKO = zeros(nDelGenes1,nDelGenes2);

for geneNo1 = 1:nDelGenes1   
    % Find gene index
    [isInModel,geneID1] = ismember(geneList1{geneNo1},model.genes);
    if (~differentSetsFlag)
        %grRateKO(geneNo1,geneNo1) = singleRate(geneNo1);
        initID = geneNo1+1;
    else
        initID = 1;
    end
    parfor geneNo2 = initID:nDelGenes2
        [isInModel,geneID2] = ismember(geneList2{geneNo2},model.genes);
        % Find rxns associated with this gene
        rxnInd = find(any(model.rxnGeneMat(:,[geneID1 geneID2]),2));
        rxnInd1 = find(any(model.rxnGeneMat(:,[geneID1]),2));
        rxnInd2 = find(any(model.rxnGeneMat(:,[geneID2]),2));
        if (~isempty(rxnInd))
            % Initialize the state of all genes
            x = true(nGenes,1);
            x(geneID1) = false;
            x(geneID2) = false;
            modelTmp = model;
            % Figure out if any of the reaction states is changed
			% Due to lack of transcriptional information and the fact that we are using KD instead of KO
			% We will assume each isozyme is responsible for 100% of the flux through the reaction
            % Mar 16,2011: Need to alter to still handle KO cases though: should be fixed now.
            % Use FBA/MOMA/lMOMA to calculate deletion strain growth rate
            if mv1 > 0
              for ri = 1:length(rxnInd1)
		if abs(WTflux(rxnInd1(ri))) > 0
		  if WTflux(rxnInd1(ri)) > 0
		     modelTmp.ub(rxnInd1(ri)) = WTflux(rxnInd1(ri))*mv1;
  		  else
		     modelTmp.lb(rxnInd1(ri)) = WTflux(rxnInd1(ri))*mv1;
		  end
                end
              end
            end
            if mv2 > 0
              for ri = 1:length(rxnInd2)
                if abs(WTflux(rxnInd2(ri))) > 0
		  if WTflux(rxnInd2(ri)) > 0
		    modelTmp.ub(rxnInd2(ri)) = WTflux(rxnInd2(ri))*mv2;
		  else
		    modelTmp.lb(rxnInd2(ri)) = WTflux(rxnInd2(ri))*mv2;
		  end
                end
              end
            end
            % Do KOs now, in the case of isozymes this will overwrite weaker bounds.
            if mv1 == 0
	      constrainRxn = false(length(rxnInd1),1);
              for ri = 1:length(rxnInd1)
                  if (~eval(model.rules{rxnInd1(ri)}))
                    constrainRxn(ri) = true;
                  else
                    if abs(WTflux(rxnInd1(ri))) > 0 
		      if WTflux(rxnInd1(ri)) > 0
	                modelTmp.ub(rxnInd1(ri)) = WTflux(rxnInd1(ri))*dlvl;
         	      else
		        modelTmp.lb(rxnInd1(ri)) = WTflux(rxnInd1(ri))*dlvl;
                      end
                    end
                  end
              end
	      if (any(constrainRxn))
                constrRxnInd = rxnInd1(constrainRxn);
                modelTmp.lb(constrRxnInd) = 0;
                modelTmp.ub(constrRxnInd) = 0;
              end
            end
            if mv2 == 0
	      constrainRxn = false(length(rxnInd2),1);
              for ri = 1:length(rxnInd2)
                  if (~eval(model.rules{rxnInd2(ri)}))
                    constrainRxn(ri) = true;
                  else
                    if abs(WTflux(rxnInd2(ri))) > 0
	              if WTflux(rxnInd2(ri)) > 0
		        modelTmp.ub(rxnInd2(ri)) = WTflux(rxnInd2(ri))*dlvl;
		      else
		        modelTmp.lb(rxnInd2(ri)) = WTflux(rxnInd2(ri))*dlvl;
                      end
                    end
                  end
              end
	      if (any(constrainRxn))
                constrRxnInd = rxnInd2(constrainRxn);
                modelTmp.lb(constrRxnInd) = 0;
                modelTmp.ub(constrRxnInd) = 0;
              end
            end

	    %Consider case where genes are isozymes (or in same complex)
            if mv1 > 0 && mv2 > 0
	    rxnInt = intersect(rxnInd1, rxnInd2);
			   for ri = 1:length(rxnInt)
                            if abs(WTflux(rxnInt(ri))) > 0
			      if WTflux(rxnInt(ri)) > 0
				      modelTmp.ub(rxnInt(ri)) = WTflux(rxnInt(ri))*(mv2*mv1);  %Change to multiplicative effect instead of average
		   	      else
				      modelTmp.lb(rxnInt(ri)) = WTflux(rxnInt(ri))*(mv2*mv1);  %Change to multiplicative effect instead of average
			      end
                            end
			   end
            end
            solKO = optimizeCbModel(modelTmp,'max');
            solKO.f;
            if (solKO.stat > 0)
              grRateKO(geneNo1,geneNo2) = solKO.f;
              %grRateKO(geneNo2,geneNo1) = solKO.f; 
            else
              grRateKO(geneNo1,geneNo2) = 0;
              %grRateKO(geneNo2,geneNo1) = 0;
            end
        end
        %if (differentSetsFlag)
        %    grRateKO(geneNo2,geneNo1) = grRateKO(geneNo1,geneNo2);
        %end
    end
end

grRateKOTens(:,:,mval1,mval2) = grRateKO;

end
end
