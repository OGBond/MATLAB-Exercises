function sol = FDDirichlet(f,N,a,b,A,B)

    % This solves the differential equation u''(x) = f(x) subject to
    % boundary conditions u(a) = A, u(b) = B. Here, u_i means u(x(i)) where
    % x(i) = a + i*dx and dx = (b-a)/N.

    fVals = zeros(N+1,1);   % storage space for values of f(x)
    
    dx = (b-a)/N;           % setting step size of x
    
    xVals = a:dx:b;          % creates set of x-values over which the ODE 
                            % is to be solved
    
    bVals = zeros(N,1);     % creates storage space for boundary values
    bVals(1) = A;           % sets u(x(0)) = A
    bVals(N) = B;           % sets u(x(N)) = B
    
    for i = 1:N + 1
        fVals(i) = f(a + i*dx);
                            % evaluates values of f(x) over [a,b]
    end
    
    FirstRowOfD2 = [-2, 1, zeros(1, N - 1)]; 
    % creates a typical row to be put into the matrix so that central
    % difference approximations to u''(x) are calculated
    
    FirstRowOfD2 = sparse(FirstRowOfD2);
    % creates a sparse matrix to include only the non-zero entries, and
    % their coordinates in the matrix D2
    
    D2 = toeplitz(FirstRowOfD2);
    % Creates a Toeplitz matrix using this vector, i.e. one where entries 
    % on diagonals are constants, ie. where D2(i,j) = D2(i+1,j+1) for all i,j
    
    D2(1,1) = 1;      % Changes first row to create an equation u_0 = A
    D2(N+1,N+1) = 1;  % Changes final row to create an equation u_N = B
    D2(N+1,N) = 1;
    D2 = D2/dx^2;      % Divides by h^2 (as in central difference formula)
    
    % keyboard
    
    sol = D2\fVals;   % Solves the matrix equation D2*u = f for u
    plot(xVals,sol)   % Creates a plot of the solution
    
    xlabel('$x$','Interpreter','latex')
    ylabel('$u(x)$','Interpreter','latex')
    
    ax = gca;
    ax.FontSize = 14;
    
    title('Solution to $u^{\prime\prime}(x)=f(x)$ on $x\in [a,b]$','Interpreter','LaTeX')