function bounds = bounds(dim, min, max)
    bounds=ones(dim,2);
    bounds(:,1)=bounds(:,1)*min;
    bounds(:,2)=bounds(:,2)*max;
end