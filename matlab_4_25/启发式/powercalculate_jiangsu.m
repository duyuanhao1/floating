function power = powercalculate_weibull(wind,address,prob)
%powercalculate �����糡���ʣ�����Ϊ���λ������򡢷��٣�
% ��ʽΪwind(v1^3 ang1 p1;v2^3 ang2 p2;.....)��address(x1 y1;x2 y2...��
% ���Ϊ��糡����power
DD=Wake_effect_v4(address,wind);
WW=DD.^2;  %vvΪ��Ӧ���ٶȱ���
I=size(address,1); %IΪ��Ҫѡȡ�ĵ�ַ����ɢ����,��ʽΪaddress(x1 y1;x2 y2...��
D=size(wind,1);%DΪ�����ɢ����,��ʽΪwind(v1 ang1 p1;v2 ang2 p2;.....�����Ƕ�Ϊ����������˳ʱ�����ĽǶȴ�С
power_rate=zeros(size(address,1),size(wind,1));
power_temp=0;Nt=40;p_rate=5000;
M=20000;diameter=58*2;constraintviolation=0;
for j=2:I
       for k=1:j-1
          if(norm(address(j,:)-address(k,:))<5*diameter) %��ȫ����Ϊ5D    
              constraintviolation=constraintviolation+1;
              break %������һ��ѭ��
          end
       end
%       if(constraintviolation~=0) %��ȫ����Ϊ5D    
%           break  %�����ڶ���ѭ��
%       end
end
for dd=1:D
    for jj= 1:I
        power_rate(jj,dd)=(1-sqrt(sum(WW(:,jj,dd))))^3;  %��v^3���㹦��;
    end
    power_temp=power_temp+2.8935*wind(dd,1)*sum(power_rate(:,dd));
end
power=power_temp-M*constraintviolation+Nt*p_rate*sum(prob(:,1));
end

