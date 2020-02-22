function G = Gz(f,d,is_Gz)
    if is_Gz
        G = permute(f, [4 2 3 1]);
        for i = 1:length(d)
            G(:,:,:,i) = G(:,:,:,i) * d(i);
        end
        G = permute(G, [4 2 3 1]);
    else
        G = f;
    end
end