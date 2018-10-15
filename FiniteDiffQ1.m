% Finite Difference Problems

% Question 1

f = @(x) 10*sin(20*x) + cos(x^5);

FDDirichlet(f,100,0,1,0,0.1)