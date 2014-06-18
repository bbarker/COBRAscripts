%growPerc = 20;
%
% if undefined, instead of using growPerc*WT/100 as the
% growth rate, use a 'maximum' growthrate for that
% condition.

%
mutLevel = 0.5;
%

mmEthNone = mmEth;
mmEthNone.lb(536)=0;
mmNone = mm;
mmNone.lb(550)=0;
mods = {mmNone, mmEthNone, mm};

lowerBoundMax = 1000;
mmsolWT = optimizeCbModel(mm, 'max', 'one');


%Uptake abbreviates:
abbr = {'Eth', 'Glyc', 'SAdeMeth', 'Acetald', 'Acetate', 'Trehalose', ...
	'Xanthosine', 'Glutamine', 'Glutamate', 'AdeBisph','Alanine', 'Adenosine', ...
	'Allantoin', 'Arginine', 'LowGlucose', 'Phos', 'High Glucose'};

modID = 1 + [1, 1, 0, 1, 1, 0, ...
	 0, 1, 1, 0, 1, 0, ...
	 1, 1, 1, 2];

uptRxn = [536, 555, 507, 499, 498, 632, ...
	  641, 551, 552, 601, 504, 502, ...
	  505, 510, 550, 607];


SDMat = zeros(length(mm.genes), length(uptRxn)+1);
FluxMat = zeros(length(mm.rxns), length(uptRxn)+1);

for i = 1:(length(uptRxn) + 1)
  uptBnd = nan;
  if i <= length(uptRxn)
    rxn = uptRxn(i)
    mID = modID(i)
    met = abbr(i)
    mmTMP = mods{mID};
    %disp([met rxn]);
    if exist('growPerc', 'var')
      uptBnd = uptakeByGrowth(mmTMP, [growPerc], rxn, mmsolWT.f);
      mmTMP.lb(rxn) = uptBnd;
      soltmp = optimizeCbModel(mmTMP, 'max', 'one');
    else
      mmTMP.lb(rxn) = -lowerBoundMax;
      soltmp = optimizeCbModel(mmTMP, 'max', 'one');
      uptBnd = soltmp.x(rxn);
      %if a high-energy molecule, constrain appropriately:
      if soltmp.f > 1.1 * mmsolWT.f
        uptBnd = uptakeByGrowth(mmTMP, [100], rxn, mmsolWT.f);
        mmTMP.lb(rxn) = uptBnd;
        soltmp = optimizeCbModel(mmTMP, 'max', 'one');
      end
    end
  else
    mmTMP = mm;
    soltmp = mmsolWT;
    rxn = 607;
    uptBnd = mmsolWT.x(rxn);
  end
  if length(soltmp.x) > 1
    savename = ['grSC' abbr{i} num2str(mutLevel) 'x' num2str(mutLevel)];
    variedUptakeEpiSim(mmTMP,'FBA', savename, rxn, [uptBnd], mutLevel);
  end
end


