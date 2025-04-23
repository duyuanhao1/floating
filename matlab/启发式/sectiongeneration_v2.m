function layout_section = sectiongeneration_v2(address,site,gridlength,gridwidth)
%   sectiongeneration 生成可用于第二阶段规划的布局区间
%   输入初始地址address,风机布局二进制变量site，网格长度与网格宽度girdlength, girdwidth
turbine_address=[];layout_section=[];
for n=1:length(site)
    if((0.999<=site(n)&&site(n)<=1.001))
        turbine_address=[turbine_address;address(n,1),address(n,2)];
    end
end
Nt=length(turbine_address);
addresslineform=reshape(turbine_address,1,[]);
layout_section=[addresslineform(1:Nt)-0.5*gridlength,addresslineform(Nt+1:2*Nt)-0.5*gridlength;addresslineform(1:Nt)+1.5*gridlength,addresslineform(Nt+1:2*Nt)+gridwidth];
for n=1:numel(layout_section)
    if(layout_section(n)<0)
        layout_section(n)=0;
    end
    if(layout_section(n)>2000)
        layout_section(n)=2000;
    end
end
end

