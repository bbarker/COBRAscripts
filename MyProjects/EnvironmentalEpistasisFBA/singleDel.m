growPerc = 20;
mmEthNone = mmEth;
mmEthNone.lb(536)=0;
mmNone = mm;
mmNone.lb(550)=0;
mods = {mmNone, mmEthNone, mm};

mmsolWT = optimizeCbModel(mm,'max','one');

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

cell2csv('SingleDeletions.csv',abbr,',',2000);

for i = 1:(length(uptRxn) + 1)
  if i <= length(uptRxn)
    rxn = uptRxn(i)
    mID = modID(i)
    met = abbr(i)
    mmTMP = mods{mID};
    %disp([met rxn]);
    uptBnd = uptakeByGrowth(mmTMP,[growPerc],rxn,mmsolWT.f);
    mmTMP.lb(rxn) = uptBnd;
    soltmp = optimizeCbModel(mmTMP,'max','one');
  else
    mmTMP = mm;
    soltmp = mmsolWT;
  end
  if length(soltmp.x) > 1
    sdFlux = singleGeneMutation8_flux(mmTMP,'FBA',soltmp.x,0,[0]);
    FluxMat(:,i) = columnVector(soltmp.x);
    SDMat(:,i) = columnVector(sdFlux(:,1577));
  end
end

%columnVector each of these and and append to the abbr header above
dlmwrite('SingleDeletions.csv', SDMat, '-append')

