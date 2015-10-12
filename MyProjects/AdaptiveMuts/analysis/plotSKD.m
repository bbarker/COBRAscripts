function retval = plotSKD(grKD, relWTFit)

relFit = grKD/relWTFit;
[ngen nmut] = size(grKD);
nbins = 22;
binvec = [0:(nbins-1)]/10;
Ybar = zeros(nbins,nmut);

%yticks = cell(22,1);
%for i = 1:nbins
%  yticks{i} = num2str(binvec(i));
%end

for i = 1:nmut
  Ybar(:,i) = hist(relFit(:,i),binvec);
end

Ybar = log(Ybar);
bar3(Ybar);
set(gca,'XTickLabel',{'0.1','0.3','0.5','0.7','0.9','1.1','1.4','1.6','1.8','2.0'});
%set(gca,'YTickLabel',yticks);
yticks = [0 2 4 6 8 10 12 14 16 18 20 22 24 25];
yTickLabels = (yticks-1)/10;
set(gca,'YTick',yticks)
set(gca,'YTickLabel',yTickLabels)
xlabel('Fraction of WT flux in mutant gene')
ylabel('Relative Fitness')
zlabel('log(# of genes)')

