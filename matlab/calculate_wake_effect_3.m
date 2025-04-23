function DD = calculate_wake_effect_3(address, wind, prob_matrix_1,prob_matrix_2)
    % ����������
    num_big_grids = size(address, 1);
    % ÿ����������С��������
    num_small_grids_per_big = 3*3;
    column=sqrt(num_small_grids_per_big);
    % �������
    num_wind_conditions = size(wind, 1);

    % ��ʼ�� DD ����
    DD = zeros(num_big_grids, num_big_grids, num_wind_conditions);

    % β����ɢϵ��
    alpha = 0.1;
    % ����뾶
    rd = 20;
    % �������ϵ��
    CT = 0.88;

    % ������߳�
    big_grid_length = 150;
    % С����߳�
    small_grid_length = big_grid_length / column;

    % ����ÿ�ַ��
    for k = 1:num_wind_conditions
        % �������в���Ӱ��Ĵ����� A
        for B = 1:num_big_grids
            % ��ȡ������ A ��λ��
            xB = address(B, 1);
            yB = address(B, 2);

            % ����������Ӱ��Ĵ����� B
            for A = 1:num_big_grids
                % ��� A �� B ��ͬһ������������������Ӱ��Ϊ 0
                if A == B
                    DD(A, B, k) = 0;
                    continue;
                end

                % ��ȡ������ B ��λ��
                xA = address(A, 1);
                yA = address(A, 2);

                % ��Ȩ��ͳ�ʼ��
                weighted_sum = 0;

                % ���������� A �ڵ�С����
                for i = 1:num_small_grids_per_big
                    % ����С���� i �ڴ����� A �е�����ƫ��
                    [offset_xA, offset_yA] = get_small_grid_offset(i, small_grid_length, column);
                    % С���� i �ڴ����� A �е�ʵ������
                    small_xA = xA + offset_xA;
                    small_yA = yA + offset_yA;

                    % ��ȡС���� i �ķ��ø���
                    prob_Ai = prob_matrix_2(i);

                    % ���������� B �ڵ�С����
                    for j = 1:num_small_grids_per_big
                        % ����С���� j �ڴ����� B �е�����ƫ��
                        [offset_xB, offset_yB] = get_small_grid_offset(j, small_grid_length, column);
                        % С���� j �ڴ����� B �е�ʵ������
                        small_xB = xB + offset_xB;
                        small_yB = yB + offset_yB;

                        % ��ȡС���� j �ķ��ø���
                        prob_Bj = prob_matrix_1(j);

                        % ����β����ʧ
                        result = wake_loss([small_xA, small_yA], [small_xB, small_yB], wind(k, :), alpha, rd, CT);

                        % �ۼӼ�Ȩβ����ʧ
                        weighted_sum = weighted_sum + prob_Ai * prob_Bj * result;
                    end
                end

                % ����Ȩ��ͽ������ DD ����
                DD(B, A, k) = weighted_sum;
                DD = abs(real(DD));
            end
        end
    end
end

function [offset_x, offset_y] = get_small_grid_offset(small_grid_index, small_grid_length, column)
    % ����С�����ڴ������ڵ��к�������
    small_row = mod(small_grid_index - 1, column) + 1;
    small_col = floor((small_grid_index - 1) / column) + 1;
    % ����С���������ƫ��
    offset_y = (small_row - 1) * small_grid_length;
    offset_x = (small_col - 1) * small_grid_length;
end

function result = wake_loss(address_i, address_j, wind_k, alpha, rd, CT)
    % �����������ֵ��ڷ����ϵľ���
    d = (address_i(2)-address_j(2))*cos(wind_k(2))+(address_i(1)-address_j(1))*sin(wind_k(2)); 
    % �����������ֵ��ڴ�ֱ�ڷ������ϵľ���
    d_v = abs((address_i(2)-address_j(2))*sin(wind_k(2))-(address_i(1)-address_j(1))*cos(wind_k(2))); 

    if(d <= 0)  % dС�ڵ���0��β��������Ӱ��
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
    