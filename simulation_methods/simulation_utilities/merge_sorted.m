function [merged, aux] = merge_sorted(v1,v2,options)
%MERGE_SORTED merges a pair of sorted vectors, and also constructs an
%auxilliary vector if there is some accompanying calculation described by
%auxfun.
%
%last updated 09/22/26 by Adam Petrucci
arguments (Input)
    v1       % first sorted vector
    v2       % second sorted vector
end
arguments (Input)
    options.AuxFun = @(z1,z2,h1,h2) 0  % the auxiliary function
end
arguments (Output)
    merged     % the merged vector
    aux        % the auxilliary calculations
end

    auxfun = options.AuxFun;

    n1 = length(v1);
    n2 = length(v2);

    merged = zeros(1,n1+n2);
    aux = zeros(1,n1+n2);

    i1 = 1;  % counter for v1
    i2 = 1;  % counter for v2
    k  = 0;  % length of merged vector

    while i1 <= n1 || i2 <= n2

        k = k+1;  % update vector length

        if i1 > n1   % only v2 left

            merged(k) = v2(i2);
            aux(k) = auxfun(i1-1,i2,false,true);
            i2 = i2 + 1;

        elseif i2 > n2  % only v1 left

            merged(k) = v1(i1);
            aux(k) = auxfun(i1,i2-1,true,false);
            i1 = i1 + 1;

        elseif v1(i1) < v2(i2)   % append next v1 element

            merged(k) = v1(i1);
            aux(k) = auxfun(i1,i2-1,true,false);
            i1 = i1 + 1;

        elseif v2(i2) < v1(i1)   % append next v2 element

            merged(k) = v2(i2);
            aux(k) = auxfun(i1-1,i2,false,true);
            i2 = i2 + 1;

        else  % v1 == v2, append once but increment both

            merged(k) = v1(i1);
            aux(k) = auxfun(i1,i2,true,true);
            i1 = i1 + 1;
            i2 = i2 + 1;

        end

    end

    % adjust size (in case there were doubles;
    merged = merged(1:k);
    aux = aux(1:k);

end