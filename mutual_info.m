function I = mutual_info(x,y,nbins)

    xb = discretize(x,nbins);
    yb = discretize(y,nbins);

    Pxy = accumarray([xb yb], 1, [nbins nbins]) / numel(x);
    Px = sum(Pxy,2);
    Py = sum(Pxy,1);

    I = 0;
    for a = 1:nbins
        for b = 1:nbins
            if Pxy(a,b) > 0
                I = I + Pxy(a,b) * log( Pxy(a,b) / (Px(a) * Py(b)) );
            end
        end
    end

end
