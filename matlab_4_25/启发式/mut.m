function SelCh = mut(SelCh, pm)
    [Nind, Lind] = size(SelCh);
    % 遍历每个个体的每个基因位
    for i = 1:Nind
        for j = 1:Lind
            % 以变异概率 pm 决定是否进行变异
            if rand < pm
                % 变异操作：取反
                SelCh(i, j) = ~SelCh(i, j);
            end
        end
    end
end