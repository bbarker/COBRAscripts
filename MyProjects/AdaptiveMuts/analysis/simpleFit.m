function paramVec = simpleFit(xdata, ydata, x0)

  %Here are functions to try in the fit%
  %pick one and use it below%

  function F = expCon(x, xdata)
    %x(1) == a; x(2) == b; x(3) == c
    F = x(1)*(xdata^x(2)) + x(3);
  end

  %change the first argument below to specify desired 
  %function, from the list of functions above (e.g. expCon)
  %or an external function.
  paramVec = lsqcurvefit(@expCon,x0,xdata,ydata);

end