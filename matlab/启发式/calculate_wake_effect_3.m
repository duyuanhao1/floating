function DD = calculate_wake_effect_3(address, wind, prob_matrix_1,prob_matrix_2)
    % 大网格数量
    num_big_grids = size(address, 1);
    % 每个大网格内小网格数量
    num_small_grids_per_big = 3*3;
    column=sqrt(num_small_grids_per_big);
    % 风况数量
    num_wind_conditions = size(wind, 1);

    % 初始化 DD 矩阵
    DD = zeros(num_big_grids, num_big_grids, num_wind_conditions);

    % 尾流扩散系数
    alpha = 0.1;
    % 风机半径
    rd = 20;
    % 风机推力系数
    CT = 0.88;

    % 大网格边长
    big_grid_length = 150;
    % 小网格边长
    small_grid_length = big_grid_length / column;

    % 遍历每种风况
    for k = 1:num_wind_conditions
        % 遍历所有产生影响的大网格 A
        for B = 1:num_big_grids
            % 获取大网格 A 的位置
            xB = address(B, 1);
            yB = address(B, 2);

            % 遍历所有受影响的大网格 B
            for A = 1:num_big_grids
                % 如果 A 和 B 是同一个大网格，自身对自身的影响为 0
                if A == B
                    DD(A, B, k) = 0;
                    continue;
                end

                % 获取大网格 B 的位置
                xA = address(A, 1);
                yA = address(A, 2);

                % 加权求和初始化
                weighted_sum = 0;

                % 遍历大网格 A 内的小网格
                for i = 1:num_small_grids_per_big
                    % 计算小网格 i 在大网格 A 中的坐标偏移
                    [offset_xA, offset_yA] = get_small_grid_offset(i, small_grid_length, column);
                    % 小网格 i 在大网格 A 中的实际坐标
                    small_xA = xA + offset_xA;
                    small_yA = yA + offset_yA;

                    % 获取小网格 i 的放置概率
                    prob_Ai = prob_matrix_2(i);

                    % 遍历大网格 B 内的小网格
                    for j = 1:num_small_grids_per_big
                        % 计算小网格 j 在大网格 B 中的坐标偏移
                        [offset_xB, offset_yB] = get_small_grid_offset(j, small_grid_length, column);
                        % 小网格 j 在大网格 B 中的实际坐标
                        small_xB = xB + offset_xB;
                        small_yB = yB + offset_yB;

                        % 获取小网格 j 的放置概率
                        prob_Bj = prob_matrix_1(j);

                        % 计算尾流损失
                        result = wake_loss([small_xA, small_yA], [small_xB, small_yB], wind(k, :), alpha, rd, CT);

                        % 累加加权尾流损失
                        weighted_sum = weighted_sum + prob_Ai * prob_Bj * result;
                    end
                end

                % 将加权求和结果存入 DD 矩阵
                DD(B, A, k) = weighted_sum;
                DD = abs(real(DD));
            end
        end
    end
end

function [offset_x, offset_y] = get_small_grid_offset(small_grid_index, small_grid_length, column)
    % 计算小网格在大网格内的行和列索引
    small_row = mod(small_grid_index - 1, column) + 1;
    small_col = floor((small_grid_index - 1) / column) + 1;
    % 计算小网格的坐标偏移
    offset_y = (small_row - 1) * small_grid_length;
    offset_x = (small_col - 1) * small_grid_length;
end

function result = wake_loss(address_i, address_j, wind_k, alpha, rd, CT)
    % 计算两个积分点在风向上的距离
    d = (address_i(2)-address_j(2))*cos(wind_k(2))+(address_i(1)-address_j(1))*sin(wind_k(2)); 
    % 计算两个积分点在垂直于风向方向上的距离
    d_v = abs((address_i(2)-address_j(2))*sin(wind_k(2))-(address_i(1)-address_j(1))*cos(wind_k(2))); 

    if(d <= 0)  % d小于等于0，尾流对其无影响
        result = 0;
    else     % d > 0
        if(d_v <= d*alpha)
            result = (1-sqrt(1-CT))*((rd/(rd+alpha*d))^2);
        elseif(d*alpha < d_v && d_v <= rd + rd + d*alpha)  
            A_shadow=acos(((d*alpha+rd)^2+d_v^2-rd^2)/(2*(d*alpha+rd)*d_v))*(d*alpha+rd)^2+acos((-(d*alpha+rd)^2+d_v^2+rd^2)/(2*(rd)*d_v))*(rd)^2-sin(acos(((d*alpha+rd)^2+d_v^2-rd^2)/(2*(d*alpha+rd)*d_v)))*(d*alpha+rd)*d_v;
            ita=A_shadow/(pi*rd^2);
            result=(1-sqrt(1-CT))*((rd/(rd+alpha*d))^2)*ita;   
        elseif(d_v > rd + rd + d*alpha)  
            result = 0;
        end
    end
end
    