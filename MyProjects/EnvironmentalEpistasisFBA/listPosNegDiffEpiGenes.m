function listPosNegDiffEpiGenes()

%
% Configure the below lines to suit your needs
%
diffEpiCutoff = 0.01;
WTFile    = '/home/brandon/FBA/models/Analysis/uptakeGradients/WT.txt';
condFile  = '/home/brandon/FBA/models/Analysis/uptakeGradients/SmallScale20/epigrSCEth0p5x0p5-1.9159';
genesFile = '/home/brandon/FBA/models/mmgenes.csv';
%
% End of config
%

genes = importdata(genesFile);

WT_epi     = dlmread(WTFile); 
cond_epi   = dlmread(condFile);

dimsEq = length(genes)
dimsEq = union(dimsEq, size(WT_epi))
dimsEq = union(dimsEq, size(cond_epi))
if length(dimsEq) ~= 1
%
disp(['Error: some dimensions are not equal for' ... 
  ' what should be the number of genes.']);
return;
%
end

epiDiff    = cond_epi - WT_epi;
posDEgenesMap = containers.Map;
negDEgenesMap = containers.Map;

for i = 1:dimsEq
for j = 1:dimsEq
%
if i ~= j
if epiDiff(i, j) >= diffEpiCutoff
  posDEgenesMap(genes{i}) = 1;
  posDEgenesMap(genes{j}) = 1;
elseif epiDiff(i, j) <= -diffEpiCutoff
  negDEgenesMap(genes{i}) = 1;
  negDEgenesMap(genes{j}) = 1;
else
  1==1;
end
end % if i ~= j
%
end % for j
end % for i


posDEgenes = keys(posDEgenesMap);
posDEgenes = posDEgenes';
negDEgenes = keys(negDEgenesMap);
negDEgenes = negDEgenes';
szpos = size(posDEgenes)
szneg = size(negDEgenes)

% Unsafe writing ...
cell2csv('posDiffEpiGenes.csv', posDEgenes, ',');
cell2csv('negDiffEpiGenes.csv', negDEgenes, ',');
