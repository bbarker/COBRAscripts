function uptakeVec = uptakeByGrowth(model,growthPerc,limitNutrient,maxfit,otype)
%uptakeByGrowth Returns the uptake rates of a nutrient that constrain the model to the given fitnesses, relative to WT
%INPUT
% model         COBRA model structure including gene-reaction associations
% growthPerc    Vector of growth percentages relative to WT=100%
% limitNutrient the position of the uptake rate of a given nutrient in the flux vector; could also be a generic reaction.
%
%OUTPUTS
% uptakeVec       Vector of uptakeRates

solWT = optimizeCbModel(model,'max','one');
if nargin < 4
  maxfit = solWT.f;
end

if nargin < 5
  otype = 'max';
end

%assume biomass or single objective here,
%since we would need to directly call a CobraSolver to 
%add constraints on generic prior objectives
bpos = find(model.c);
uptakeVec = zeros(length(limitNutrient),length(growthPerc));



for i = 1:length(growthPerc)
  mtmp = model;
  mtmp.c(bpos)=0;
  mtmp.lb(bpos)=maxfit*growthPerc(i)/100;
  if length(limitNutrient) > 1
    for j = 1:length(limitNutrient)
      mtmp.c(limitNutrient(j)) = 1/solWT.x(limitNutrient(j)); %why is this 1/limitNutrient? Shouldn't really matter
      mtmp.lb(limitNutrient(j)) = -1000;
    end
  else
    mtmp.c(limitNutrient)=1;
    mtmp.lb(limitNutrient) = -1000;
  end
  %depending on exchange reaction format in model, may need to minimize or maximize.
  solt = optimizeCbModel(mtmp,otype);
  if length(solt.x) > 1
    if length(limitNutrient) > 1
      for j = 1:length(limitNutrient)
        uptakeVec(j,i) = solt.x(limitNutrient(j));
      end
    else
      tmp = solt.x(limitNutrient);
      uptakeVec(i) = solt.x(limitNutrient);
    end
  end
end
