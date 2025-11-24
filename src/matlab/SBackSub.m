function x = SBackSub(U,d)
    if iscell(U) == false
        x = BackSub(U,d);
        return
    end

    n = size(U,1);
    x = zeros(n,1);

    x(n) = d(n)/U{n}(n);

    for i = n-1:-1:1
        pivot_recip = 1/U{i}(i);

        for j = 1:numEntries(U{i})
            all_keys = keys(U{i});
            jth_key = all_keys(j);

            x(i) = x(i)+U{i}(jth_key)*x(jth_key);
        end
        
        x(i) = (d(i)-x(i))*pivot_recip;
    end
end