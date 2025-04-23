function SelCh = select(method, Chrom, FitnV, GGAP)
    if strcmp(method, 'sus')
        [Nind, ~] = size(Chrom);
        NSel = round(GGAP * Nind); % 要选择的个体数量
        cumFitnV = cumsum(FitnV); % 计算累积适应度值
        step = cumFitnV(end) / NSel; % 步长
        start = rand * step; % 随机起始点
        pointers = start:step:start + (NSel - 1) * step; % 选择指针
        SelCh = zeros(NSel, size(Chrom, 2));
        i = 1;
        j = 1;
        while i <= NSel
            if pointers(i) <= cumFitnV(j)
                SelCh(i, :) = Chrom(j, :);
                i = i + 1;
            else
                j = j + 1;
            end
        end
    else
        error('Unsupported selection method.');
    end
end