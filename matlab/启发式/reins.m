function [Chrom, ObjV] = reins(Chrom, SelCh, Nkeep, Nmult, ObjV, ObjVSel)
    [Nind_chrom, ~] = size(Chrom);
    [Nind_selch, ~] = size(SelCh);
    
    % 按目标函数值对父代种群进行排序，获取索引
    [~, sort_indices_chrom] = sort(ObjV, 'descend');
    % 按目标函数值对子代种群进行排序，获取索引
    [~, sort_indices_selch] = sort(ObjVSel, 'descend');
    
    % 保留父代中适应度最好的 Nkeep 个个体
    Chrom_best = Chrom(sort_indices_chrom(1:Nkeep), :);
    ObjV_best = ObjV(sort_indices_chrom(1:Nkeep));
    
    % 从子代中选择 Nind_chrom - Nkeep 个个体插入
    num_to_select = Nind_chrom - Nkeep;
    if Nmult > 1
        % 如果 Nmult > 1，子代个体可以重复选择
        sel_indices = randi(Nind_selch, num_to_select, 1);
    else
        % 如果 Nmult = 1，子代个体不重复选择
        sel_indices = sort(randi(Nind_selch, num_to_select, 1));
    end
    Chrom_sel = SelCh(sel_indices, :);
    ObjV_sel = ObjVSel(sel_indices);
    
    % 组合新的种群和目标函数值
    Chrom = [Chrom_best; Chrom_sel];
    ObjV = [ObjV_best; ObjV_sel];
end