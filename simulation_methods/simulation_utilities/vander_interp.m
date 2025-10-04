function return_data = vander_interp(x,b)
% I should try to figure out the 3-tensor formalism for this, it would
% streamline the code

    n = length(x);
    m = length(b)/n;

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

    %for k = 1:n
    %    x_(m*(k-1)+1:m*k) = ones([1,m])*x(k);
    %end
    %x = x_;

    %P = (1:n*m)-1;
    %XP = fliplr(x.^P);

    V = zeros(n*m);
    V_ = zeros([m,n*m]);
    for i = 0:m-1
        V_(i+1,i+1:end) = arrayfun(@(x) factorial(x)/factorial(x-i),i:(n*m)-1);
    end
    V_ = fliplr(V_);

    for k = 1:n
        V((1:m)+(k-1)*m,1:n*m) = V_;
    end

    A = V.*XP;
    coeffs = A\b.';

    return_data = @(x) polyval(coeffs,x);
    
end