function fluxKOTens = singleGeneMutation8_flux(model,method,WTflux,dlvl,fitvals1,geneList1,geneList2,enzConstraint,verbFlag)

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
if (nargin < 6)
     differentSetsFlag = false;
    geneList1 = model.genes;
else
    if (isempty(geneList1))
        geneList1 = model.genes;
    end
end

if (nargin < 7)
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
if (nargin < 8)
    enzConstraint=1;
end
if (nargin < 9)
    verbFlag = false;
end

%leave out 0 for fitvals1; need to do standard KO sim for this. 
%Also, don't need 1x1 (WT)'
if nargin < 5
  fitvals1=[0.5];
end
fitvals2=[1];

nGenes = length(model.genes);
disp(length(geneList1));
disp(length(geneList1));
%solTruMax = optimizeCbModel(model,'max');
%grRateMax = solTruMax.f;
if enzConstraint ~= 1
 model = constrainEnzFluxGlobal(model, enzConstraint);  %Need to fix this later for using taxicab norm
end

solMax = optimizeCbModel(model,'max');
grRateWT = solMax.f;
if nargin < 3
   solMax = optimizeCbModel(model,'max',true,true);
   WTflux = solMax.f;
 %else
 %  disp('Using precomputed flux vector');
end


nDelGenes1 = length(geneList1);
nDelGenes2 = length(geneList2);

%grRateKOTens = ones(nDelGenes1,nDelGenes2, 5*5)*grRateWT;
%grRateKOTens = zeros(nDelGenes1, nDelGenes2, length(fitvals1), length(fitvals2));
%just assume one fitval in this specialization for parallelization
fluxKOTens = zeros(nDelGenes1, length(model.lb));

if (differentSetsFlag)
    nTotalPairs = nDelGenes1*nDelGenes2;
else
    nTotalPairs = nDelGenes1*(nDelGenes1-1)/2;
end

% Run double deletion analysis
delCounter = 0;


mv1 = fitvals1(1);
mv2 = 1;
parfor geneNo1 = 1:nDelGenes1   
    % Find gene index
    [isInModel,geneID1] = ismember(geneList1{geneNo1},model.genes);
    if (~differentSetsFlag)
        %grRateKO(geneNo1,geneNo1) = singleRate(geneNo1);
        initID = geneNo1+1;
    else
        initID = 1;
    end
    for geneNo2 = 1:1
        delCounter = delCounter + 1;
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
            constrainRxn = false(length(rxnInd),1);
            % Figure out if any of the reaction states is changed
			% (needed for isoenzymes)
			% Due to lack of transcriptional information and the fact that we are using KD instead of KO
			% We will assume each isozyme is responsible for 100% of the flux through the reaction
			% For isozyme pairs, we will add the fitness deficit.
            %for ri = 1:length(rxnInd)
            %    if (~eval(model.rules{rxnInd(ri)}))
            %        constrainRxn(ri) = true;
            %    end
            %end
            % Use FBA/MOMA/lMOMA to calculate deletion strain growth rate
            %if (any(constrainRxn))
                modelTmp = model;
              if mv1 >=0 %mv1 > 0  %uncomment if mixing mutation models
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
              if mv1 < 0 %mv1 == 0 %uncomment if mixing mutation models
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

               solKO = optimizeCbModel(modelTmp,'max','one');
               solKO.f;
               if (solKO.stat > 0)
                   fluxKOTens(geneNo1,:) = solKO.x;
               else
                   fluxKOTens(geneNo1,:) = 0*WTflux;
               end
        end
    end
end

