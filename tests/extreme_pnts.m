function [A B] = extreme_pnts(S)
    S_sort=sort(S);
    A = S_sort(1);
    B = S_sort(end);
end
