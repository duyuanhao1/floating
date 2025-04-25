function SelCh = recombin(method, SelCh, pc)
    if strcmp(method, 'xovsp')
        [Nind, Lind] = size(SelCh);
        % 确保 Nind 是偶数
        if mod(Nind, 2) ~= 0
            % 如果 Nind 是奇数，向下取整
            Nind = Nind - 1;
        end
        % 确定需要进行交叉操作的个体对
        crossPairs = rand(Nind/2, 1) < pc;
        for i = 1:Nind/2
            if crossPairs(i)
                % 选择两个个体
                ind1 = 2 * i - 1;
                ind2 = 2 * i;
                % 随机选择交叉点
                crossPoint = randi([1, Lind - 1]);
                % 进行单点交叉
                temp = SelCh(ind1, crossPoint + 1:end);
                SelCh(ind1, crossPoint + 1:end) = SelCh(ind2, crossPoint + 1:end);
                SelCh(ind2, crossPoint + 1:end) = temp;
            end
        end
    else
        error('Unsupported recombination method.');
    end
end
