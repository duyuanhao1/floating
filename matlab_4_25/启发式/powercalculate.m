function Objective = powercalculate(wind,address)
%powercalculate �����糡���ʣ�����Ϊ���λ������򡢷��٣�
% ��ʽΪwind(v1^3 ang1 p1;v2^3 ang2 p2;.....)��address(x1 y1;x2 y2...��
% ���Ϊ��糡������ʧ

% Ϊÿ�ַ��ƽ̨���岻ͬ�ĸ��ʾ���
prob_matrix_y1=[0.5/8,0.5/8,0.5/8,0.5/8,0.5,0.5/8,0.5/8,0.5/8,0.5/8]';
prob_matrix_y2=[0,0,0,0,1,0,0,0,0]';
prob_matrix_y3=[0.8/8,0.8/8,0.8/8,0.8/8,0.2,0.8/8,0.8/8,0.8/8,0.8/8]';

% ���ú�������β��Ӱ��
DD1_1 = calculate_wake_effect_3(address, wind, prob_matrix_y1,prob_matrix_y1);%y1���ƽ̨��y2���ƽ̨��β��Ӱ��
DD1_2 = calculate_wake_effect_3(address, wind, prob_matrix_y1,prob_matrix_y2);
DD1_3 = calculate_wake_effect_3(address, wind, prob_matrix_y1,prob_matrix_y3);

DD2_1 = calculate_wake_effect_3(address, wind, prob_matrix_y2,prob_matrix_y1);
DD2_2 = calculate_wake_effect_3(address, wind, prob_matrix_y2,prob_matrix_y2);
DD2_3 = calculate_wake_effect_3(address, wind, prob_matrix_y2,prob_matrix_y3);

DD3_1 = calculate_wake_effect_3(address, wind, prob_matrix_y3,prob_matrix_y1);
DD3_2 = calculate_wake_effect_3(address, wind, prob_matrix_y3,prob_matrix_y2);
DD3_3 = calculate_wake_effect_3(address, wind, prob_matrix_y3,prob_matrix_y3);

I = size(address, 1); % IΪ��Ҫѡȡ�ĵ�ַ����ɢ����,��ʽΪaddress(x1 y1;x2 y2...��
D = size(wind, 1); % DΪ�����ɢ����,��ʽΪwind(v1 ang1 p1;v2 ang2 p2;.....�����Ƕ�Ϊ����������˳ʱ�����ĽǶȴ�С

vv1_1 = 1 - DD1_1;  % vvΪ��Ӧ���ٶȱ���
vv1_2 = 1 - DD1_2;
vv1_3 = 1 - DD1_3;
vv2_1 = 1 - DD2_1;
vv2_2 = 1 - DD2_2;
vv2_3 = 1 - DD2_3;
vv3_1 = 1 - DD3_1;
vv3_2 = 1 - DD3_2;
vv3_3 = 1 - DD3_3;

loss_num1_1 = 1 - vv1_1.^3;
loss_num1_2 = 1 - vv1_2.^3;
loss_num1_3 = 1 - vv1_3.^3;
loss_num2_1 = 1 - vv2_1.^3;
loss_num2_2 = 1 - vv2_2.^3;
loss_num2_3 = 1 - vv2_3.^3;
loss_num3_1 = 1 - vv3_1.^3;
loss_num3_2 = 1 - vv3_2.^3;
loss_num3_3 = 1 - vv3_3.^3;
site= readmatrix('4_24_9101036.txt');  
y1=site(:,1);
y2=site(:,2);
y3=site(:,3);
M=50000;diameter=40;constraintviolation=0;
for j=1:I
       for k=1:j-1
          if(norm(address(j,:)-address(k,:))<=280&&y1(k) + y1(j)==2) %��ȫ����Ϊ5D    
              constraintviolation=constraintviolation+1;
              break %������һ��ѭ��
          end
          if(norm(address(j,:)-address(k,:))<=240&&y1(k) + y2(j)==2) %��ȫ����Ϊ5D    
              constraintviolation=constraintviolation+1;
              break %������һ��ѭ��
          end
          if(norm(address(j,:)-address(k,:))<=360&&y3(k) + y3(j)==2) %��ȫ����Ϊ5D    
              constraintviolation=constraintviolation+1;
              break %������һ��ѭ��
          end
          if(norm(address(j,:)-address(k,:))<=280&&y3(k) + y2(j)==2) %��ȫ����Ϊ5D    
              constraintviolation=constraintviolation+1;
              break %������һ��ѭ��
          end
          if(norm(address(j,:)-address(k,:))<=320&&y3(k) + y1(j)==2) %��ȫ����Ϊ5D    
              constraintviolation=constraintviolation+1;
              break %������һ��ѭ��
          end          
          if(norm(address(j,:)-address(k,:))<=200&&y2(k) + y2(j)==2) %��ȫ����Ϊ5D    
              constraintviolation=constraintviolation+1;
              break %������һ��ѭ��
          end
       end
end

site1=sum(site,2);
power = 0;
for dd = 1:D
    for ii = 1:I
        loss = 0;
        for jj = 1:I           
            loss = loss + ((y1(ii)) * ((y1(jj)) * loss_num1_1(ii, jj, dd) + (y2(jj)) * loss_num1_2(ii, jj, dd) + (y3(jj)) * loss_num1_3(ii, jj, dd)) + ...
                         (y2(ii)) * ((y1(jj)) * loss_num2_1(ii, jj, dd) + (y2(jj)) * loss_num2_2(ii, jj, dd) + (y3(jj)) * loss_num2_3(ii, jj, dd)) + ...
                         (y3(ii)) * ((y1(jj)) * loss_num3_1(ii, jj, dd) + (y2(jj)) * loss_num3_2(ii, jj, dd) + (y3(jj)) * loss_num3_3(ii, jj, dd)));
        end
        power = power + 2.8935 * wind(dd, 1) * wind(dd, 3) * (1 - loss) * site1(ii) - M*constraintviolation;
    end 
end
p_rate=5000;
wake_loss = 3000*0.93*25*(p_rate*(sum(y1)+sum(y2)+sum(y3))-power);
%wake_loss=-sum(p,'all');
fp_cost = sum(y1)*246800000+sum(y2)*288100000+sum(y3)*282200000;
Objective = wake_loss + fp_cost;

fprintf('����Ϊ %d\n', Objective);
end