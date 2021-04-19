function f = EvaluateIndividual(x)

fNumerator1 = exp(-x(1)^2);
fNumerator2 = sqrt(5)*(sin(x(2)*x(1)*x(2))^3);
fNumerator3 = 2*(cos(2*x(1)+3*x(2))^3);

fDenominator = 1 + x(1)^2 + x(2)^2;

f = (fNumerator1 + fNumerator2 + fNumerator3)/fDenominator;


% r = sqrt (x(1) .^ 2 + x(2) .^ 2) ;
% f = sin (r) ./ r;

end
