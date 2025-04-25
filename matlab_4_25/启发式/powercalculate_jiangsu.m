function power = powercalculate_weibull(wind,address,prob)
%powercalculate 计算风电场功率，输入为风机位置与风向、风速，
% 格式为wind(v1^3 ang1 p1;v2^3 ang2 p2;.....)，address(x1 y1;x2 y2...）
% 输出为风电场功率power
DD=Wake_effect_v4(address,wind);
WW=DD.^2;  %vv为对应的速度比例
I=size(address,1); %I为将要选取的地址的离散集合,格式为address(x1 y1;x2 y2...）
D=size(wind,1);%D为风的离散集合,格式为wind(v1 ang1 p1;v2 ang2 p2;.....），角度为与正北方向按顺时针计算的角度大小
power_rate=zeros(size(address,1),size(wind,1));
power_temp=0;Nt=40;p_rate=5000;
M=20000;diameter=58*2;constraintviolation=0;
for j=2:I
       for k=1:j-1
          if(norm(address(j,:)-address(k,:))<5*diameter) %安全距离为5D    
              constraintviolation=constraintviolation+1;
              break %跳出第一层循环
          end
       end
%       if(constraintviolation~=0) %安全距离为5D    
%           break  %跳出第二层循环
%       end
end
for dd=1:D
    for jj= 1:I
        power_rate(jj,dd)=(1-sqrt(sum(WW(:,jj,dd))))^3;  %用v^3计算功率;
    end
    power_temp=power_temp+2.8935*wind(dd,1)*sum(power_rate(:,dd));
end
power=power_temp-M*constraintviolation+Nt*p_rate*sum(prob(:,1));
end

