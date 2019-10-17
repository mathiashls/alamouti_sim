clear all
clc
close all
M = 2;
Fd = 200;
Rs = 1e4;
SNR = [0:40];

info = randint(250000,1);
info_mod = pskmod(info,M);
 
% geração dos canais Rayleigh 
h_1 = rayleighchan(1/Rs,Fd);
h_1.StoreHistory = 1;
h_2 = rayleighchan(1/Rs,Fd);
h_2.StoreHistory = 1;

info_rx_1 = filter(h_1,info_mod);
info_rx_2 = filter(h_2,info_mod);


h_1_conj = conj(h_1.PathGains);
h_2_conj = conj(h_2.PathGains);

%% SA E  MRC
for snr = SNR
    info_awgn1 = awgn(info_rx_1, snr);
    info_awgn2 = awgn(info_rx_2, snr);
    
%--------------------------------------SA    
    %Equalização
    info_ray_eq_1 = info_awgn1./h_1.PathGains;
    info_ray_eq_2 = info_awgn2./h_2.PathGains;
    
    compara_a = abs(h_1.PathGains)>abs(h_2.PathGains);
    compara_b = abs(h_2.PathGains)>abs(h_1.PathGains);
    
    res1 = compara_a.*info_ray_eq_1;
    res2 = compara_b.*info_ray_eq_2;
    
    info_ray_eq = res1+res2;
    
    info_ray_demod = pskdemod(info_ray_eq,M);
    
    [num, taxa_sa(snr+1)] = biterr(info, info_ray_demod);    

%---------------------------------------------MRC    
    rcv_1 = h_1_conj .* info_awgn1;
    rcv_2 = h_2_conj .* info_awgn2;
    
    info_rcv = rcv_1 + rcv_2;
    
    info_ray_demod = pskdemod(info_rcv,M);
    
    [num, taxa_MRC(snr+1)] = biterr(info, info_ray_demod);    
    
end

%% Alamouti

info_Impar = info_mod(1:2:end);%s1,s3,s5...
info_Par = info_mod(2:2:end);%s2,s4,s6...
 
info1_Antena0 = zeros(1, length(info_mod));
info1_Antena0(1:2:end) = info_Impar;            %Antena 0 no T=t
info1_Antena0(2:2:end) = (-1*conj(info_Par));   %Antena 0 no T=t+1
 
info2_Antena1 = zeros(1, length(info_mod));
info2_Antena1(1:2:end) = info_Par;              %Antena 1 no T=t
info2_Antena1(2:2:end) = (conj(info_Impar));    %Antena 1 no T=t+1

info1_canal_h1 = filter(h_1, info1_Antena0);
info1_canal_h2 = filter(h_2, info2_Antena1);

PG_h1 = h_1.PathGains;
PG_h2 = h_2.PathGains;

rx = info1_canal_h1+info1_canal_h2;
for snr = SNR
    %Aqui cada canal teria um awgn diferente! Piora bastante o desempenho!
    %info_awgn1 = awgn(info1_canal_h1, snr);
    % info_awgn2 = awgn(info1_canal_h2, snr);
    %  rt = info_awgn1 + info_awgn2;
  
    rt = awgn(rx, snr);
    
    s_impar = transpose(rt(1:2:end) .* transpose(conj(PG_h1(1:2:end)))) + (transpose(conj(rt(2:2:end))) .* PG_h2(2:2:end));
    s_par = transpose(rt(1:2:end) .* transpose(conj(PG_h2(1:2:end)))) - (transpose(conj(rt(2:2:end))).*PG_h1(2:2:end));

    S = zeros(1, length(info_mod));
    S(1:2:end) = s_impar;
    S(2:2:end) = s_par;
    
    info_demod = pskdemod(S,M);
    
    [num, taxa_alamouti(snr+1)] = biterr(info, transpose(info_demod));    
end
%%
figure
semilogy(SNR,taxa_sa)
hold on
semilogy(SNR,taxa_MRC)
hold on
semilogy(SNR,taxa_alamouti);
legend('AS', 'MRC','Alamouti')
title('AS (1Tx, 2Rx), MRC (1Tx, 2Rx), Alamouti (2Tx, Rx).')
xlabel('SNR')
ylabel('BER')



