function spvec = ShadowPricesHeuristic(model, flux)
spvec = zeros(1,length(model.lb));
solWT = optimizeCbModel(model,'max',0.1);
for i = 1:length(flux)
  if abs(flux(i))>1e-7
    mtmp = model;
    mtmp.lb(i)=flux(i)*0.99;
    mtmp.ub(i)=flux(i)*0.99;
    soldec = optimizeCbModel(mtmp);
    mtmp = model;
    mtmp.lb(i)=flux(i)*1.01;
    mtmp.ub(i)=flux(i)*1.01;
    solinc = optimizeCbModel(mtmp);
    sols = [solinc.f soldec.f];
    spvec(i) = 100*mean(abs(sols(abs(sols)>1e-7) - solWT.f))/solWT.f;
  end
end
