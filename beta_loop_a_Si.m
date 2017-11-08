lam1 = 1.1071;
ep = 11.6858;
A = 0.939816;
B = 8.10461e-3;
lams1 = 0.06144821;
lams2 = 0.1106997;
lams3 = 17.92656;
A1 = 1.023798;
A2 = 1.058264;
A3 = 5.280792;
lux = 299792458;
wavelength_start = 1.2;
wavelength_end = 2;
size = int16((wavelength_end - wavelength_start) / 0.01 + 1);
betas = zeros(size, 2);
count = 0;

for lam = wavelength_start : 0.01 : wavelength_end
    n_Si = 1.0 * sqrt(ep + A / lam ^ 2 + B * lam1 ^ 2 / (lam ^ 2 - lam1 ^ 2));
    %n_Si = 1.035 * sqrt(1 + 10.66842933 * lam ^ 2 / (lam ^ 2 - 0.3015116485 ^ 2) + 0.003043475 * lam ^ 2 / (lam ^ 2 - 1.13475115 ^ 2) + 1.54133408 * lam ^ 2 / (lam ^ 2 - 1104 ^ 2));
    n_SiO2 = 1.0 * sqrt(1 + (0.6961663 * lam ^ 2) / (lam ^ 2 - 0.0684043 ^ 2) + (0.4079426 * lam ^ 2) / (lam ^ 2 - 0.1162414 ^ 2) + (0.8974794 * lam ^ 2) / (lam ^ 2 - 9.896161 ^ 2));
    n_CaO = (1.83 + n_SiO2) / 2;
    count = count + 1;
    info = mphsolutioninfo(TPA3(n_Si, n_SiO2, lam));
    betas(count, 1) = lam;
    betas(count, 2) = real(1i * info.sol1.map(:,1));
    lam
end

dlmwrite('TPA_c-Si-betas_compare.txt', betas, 'precision', 13, 'delimiter', ',');

omegas = zeros(size, 1);
for a = 1:size
    omegas(a) = 2 * pi * lux / (betas(a, 1) * 1e-6);
end

ad = 2;
beta1s = zeros(size - 2 * ad, 2);
beta1s(:, 1) = betas(1 + ad, 1) : 0.01 : betas(size - ad, 1);
for a = 1: size - 2 * ad
    beta1s(a, 2) = (betas(2 * ad + a, 2) - betas(a, 2)) / (omegas(2 * ad + a) - omegas(a));
end

beta2s = zeros(size - 4 * ad, 2);
beta2s(:, 1) = betas(1 + 2 * ad, 1) : 0.01 : betas(size - 2 * ad, 1);
for a = 1: size - 4 * ad
    beta2s(a, 2) = (beta1s(2 * ad + a, 2) - beta1s(a, 2)) / (omegas(3 * ad + a) - omegas(ad + a));
end

beta3s = zeros(size - 6 * ad, 2);
beta3s(:, 1) = betas(1 + 3 * ad, 1) : 0.01 : betas(size - 3 * ad, 1);
for a = 1: size - 6 * ad
    beta3s(a, 2) = (beta2s(2 * ad + a, 2) - beta2s(a, 2)) / (omegas(4 * ad + a) - omegas(2 * ad + a));
end

beta4s = zeros(size - 8 * ad, 2);
beta4s(:, 1) = betas(1 + 4 * ad, 1) : 0.01 : betas(size - 4 * ad, 1);
for a = 1: size - 8 * ad
    beta4s(a, 2) = (beta3s(2 * ad + a, 2) - beta3s(a, 2)) / (omegas(5 * ad + a) - omegas(3 * ad + a));
end

beta5s = zeros(size - 10 * ad, 2);
beta5s(:, 1) = betas(1 + 5 * ad, 1) : 0.01 : betas(size - 5 * ad, 1);
for a = 1: size - 10 * ad
    beta5s(a, 2) = (beta4s(2 * ad + a, 2) - beta4s(a, 2)) / (omegas(6 * ad + a) - omegas(4 * ad + a));
end

beta6s = zeros(size - 12 * ad, 2);
beta6s(:, 1) = betas(1 + 6 * ad, 1) : 0.01 : betas(size - 6 * ad, 1);
for a = 1: size - 12 * ad
    beta6s(a, 2) = (beta5s(2 * ad + a, 2) - beta5s(a, 2)) / (omegas(7 * ad + a) - omegas(5 * ad + a));
end

beta7s = zeros(size - 14 * ad, 2);
beta7s(:, 1) = betas(1 + 7 * ad, 1) : 0.01 : betas(size - 7 * ad, 1);
for a = 1: size - 14 * ad
    beta7s(a, 2) = (beta6s(2 * ad + a, 2) - beta6s(a, 2)) / (omegas(8 * ad + a) - omegas(6 * ad + a));
end

beta8s = zeros(size - 16 * ad, 2);
beta8s(:, 1) = betas(1 + 8 * ad, 1) : 0.01 : betas(size - 8 * ad, 1);
for a = 1: size - 16 * ad
    beta8s(a, 2) = (beta7s(2 * ad + a, 2) - beta7s(a, 2)) / (omegas(9 * ad + a) - omegas(7 * ad + a));
end

beta9s = zeros(size - 18 * ad, 2);
beta9s(: ,1) = betas(1 + 9 * ad, 1) : 0.01 : betas(size - 9 * ad, 1);
for a = 1: size - 18 * ad
    beta9s(a, 2) = (beta8s(2 * ad + a, 2) - beta8s(a, 2)) / (omegas(10 * ad + a) - omegas(8 * ad + a));
end

beta10s = zeros(size - 20 * ad, 2);
beta10s(:, 1) = betas(1 + 10 * ad, 1) : 0.01 : betas(size - 10 * ad, 1);
for a = 1: size - 20 * ad
    beta10s(a, 2) = (beta9s(2 * ad + a, 2) - beta9s(a, 2)) / (omegas(11 * ad + a) - omegas(9 * ad + a));
end

beta11s = zeros(size - 22 * ad, 2);
beta11s(:, 1) = betas(1 + 11 * ad, 1) : 0.01 : betas(size - 11 * ad, 1);
for a = 1: size - 22 * ad
    beta11s(a, 2) = (beta10s(2 * ad + a, 2) - beta10s(a, 2)) / (omegas(12 * ad + a) - omegas(10 * ad + a));
end

beta12s = zeros(size - 24 * ad, 2);
beta12s(: ,1) = betas(1 + 12 * ad, 1) : 0.01 : betas(size - 12 * ad, 1);
for a = 1: size - 24 * ad
    beta12s(a, 2) = (beta11s(2 * ad + a, 2) - beta11s(a, 2)) / (omegas(13 * ad + a) - omegas(11 * ad + a));
end

beta13s = zeros(size - 26 * ad, 2);
beta13s(:, 1) = betas(1 + 13 * ad, 1) : 0.01 : betas(size - 13 * ad, 1);
for a = 1: size - 26 * ad
    beta13s(a, 2) = (beta12s(2 * ad + a, 2) - beta12s(a, 2)) / (omegas(14 * ad + a) - omegas(12 * ad + a));
end

figure
plot(beta2s(:, 1), beta2s(:, 2))
figure
plot(beta3s(:, 1), beta3s(:, 2))
figure
plot(beta4s(:, 1), beta4s(:, 2))
figure
plot(beta5s(:, 1), beta5s(:, 2))
figure
plot(beta6s(:, 1), beta6s(:, 2))
figure
plot(beta7s(:, 1), beta7s(:, 2))
figure
plot(beta8s(:, 1), beta8s(:, 2))
figure
plot(beta9s(:, 1), beta9s(:, 2))
figure
plot(beta10s(:, 1), beta10s(:, 2))
figure
plot(beta11s(:, 1), beta11s(:, 2))
figure
plot(beta12s(:, 1), beta12s(:, 2))
figure
plot(beta13s(:, 1), beta13s(:, 2))

dlmwrite('TPA-c-Si_beta2s_compare.txt', beta2s, 'precision', 13, 'delimiter', ',');
dlmwrite('TPA-c-Si_beta3s_compare.txt', beta3s, 'precision', 13, 'delimiter', ',');
dlmwrite('TPA-c-Si_beta4s_compare.txt', beta4s, 'precision', 13, 'delimiter', ',');
dlmwrite('TPA-c-Si_beta5s_compare.txt', beta5s, 'precision', 13, 'delimiter', ',');
dlmwrite('TPA-c-Si_beta6s_compare.txt', beta6s, 'precision', 13, 'delimiter', ',');
dlmwrite('TPA-c-Si_beta7s_compare.txt', beta7s, 'precision', 13, 'delimiter', ',');
dlmwrite('TPA-c-Si_beta8s_compare.txt', beta8s, 'precision', 13, 'delimiter', ',');
dlmwrite('TPA-c-Si_beta9s_compare.txt', beta9s, 'precision', 13, 'delimiter', ',');
dlmwrite('TPA-c-Si_beta10s_compare.txt', beta10s, 'precision', 13, 'delimiter', ',');
dlmwrite('TPA-c-Si_beta11s_compare.txt', beta11s, 'precision', 13, 'delimiter', ',');
dlmwrite('TPA-c-Si_beta12s_compare.txt', beta12s, 'precision', 13, 'delimiter', ',');
dlmwrite('TPA-c-Si_beta13s_compare.txt', beta13s, 'precision', 13, 'delimiter', ',');