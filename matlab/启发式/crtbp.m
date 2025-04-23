function pop = crtbp(NIND, NVARS)
    % NIND: 个体数目
    % NVARS: 变量的总二进制位数
    pop = randi([0, 1], NIND, NVARS);
end