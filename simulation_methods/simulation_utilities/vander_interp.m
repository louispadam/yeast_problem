function return_data = vander_interp(x,b)
%VANDER_INTERP constructs a interpolative polynomial given the desired
%derivatives at a collection of polynomials. Check 'Confluent Vandermonde
%Matrices' on Wikipedia for more information.
%
%last updated 10/08/25 by Adam Petrucci
arguments (Input)
    x       % Collection of points to match
    b       % Desired values of functions and derivatives, in form
            % [function value at x1, derivative value at x1, ...
            %  function value at xn, derivative value at xn].'
end

    % number of points to match
    n = length(x);

    % number of derivatives at each point
    m = length(b)/n;

    % Construct Vandermonde matrix, with correction for x=0
    XP = zeros(n*m);
    P = (1:n*m)-1;
    for k = 1:n
        if x(k) ~= 0
            a = ones([m,1]);
            a = a.*(x(k).^P);
            XP(m*(k-1)+1:m*k,:) = a;
        else
            a = zeros([m,n*m]);
            for i = 1:m
                a(i,i) = 1;
            end
            XP(m*(k-1)+1:m*k,:) = a;
        end
    end
    XP = fliplr(XP);

    % Construct Taylor matrix
    V = zeros(n*m);
    V_ = zeros([m,n*m]);
    for i = 0:m-1
        V_(i+1,i+1:end) = arrayfun(@(x) factorial(x)/factorial(x-i),i:(n*m)-1);
    end
    V_ = fliplr(V_);

    for k = 1:n
        V((1:m)+(k-1)*m,1:n*m) = V_;
    end

    % Solve linear system for polynomial coefficients
    A = V.*XP;
    coeffs = A\b.';

    return_data = @(x) polyval(coeffs,x);
    
end