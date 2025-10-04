
clc;
clear;
close all;

% Parameters
M = 8; % Charlie's Constellation size
alpha = [10^-2 10^-1 0.2:0.05:0.9-0.05 0.9:10^-3:1-10^-3];% Energy Splitting factor
L=length(alpha);
symbols = 100;% Num of symbols
snr_dB=35;
snr = 10^(snr_dB / 10);
Ec=1; % Charlie Avg Energy
N0 = Ec / snr; % Noise Power
N = 4; % No. of Rx Antennas at Charlie
lambda = 10^-4;% Residual SelfInterference 
sigma_AC = 2; 
iterations=1000;


%Th_P01 = zeros(1, L);
%Th_P10 =zeros(1, L);
Bob_BER=zeros(N, L);
%Charlie_BER = zeros(N, L);
Charlie_p10= zeros(N, L);
Charlie_p01 = zeros(N, L);
Omega_0=zeros(1, L); % Energy Rxd at Charlie corresponding to symbol 0
Omega_1=zeros(1, L); % Energy Rxd at Charlie corresponding to symbol 1
Threshold=zeros(1, L); % Detection Threshold at Charlie

function result=compute_BER(num_bits, var)
        result=var/num_bits;
end


for Nc = 1:N

for k = 1:L
    
    % Threshold Computation
    Omega_0(k) = N0 + lambda * alpha(k) * Ec;
    Omega_1(k) = N0 + lambda * alpha(k) * Ec + (1 - alpha(k)) * (sigma_AC^2); 
    Threshold(k) = (Nc *Omega_0(k) * Omega_1(k)) * log(Omega_0(k) / Omega_1(k)) / (Omega_0(k) - Omega_1(k));

    % P_ik Theretical Probability of decoding Alice symbol i as k
    %Th_P01(k) = (gammainc(Threshold(k)/Omega_0(k),Nc,"upper"));
    %Th_P10(k) =(gammainc(Threshold(k)/Omega_1(k),Nc,"lower")); 


N0b=N0;
N1b=N0+(1-alpha(k))*Ec;
Ec_fCB = alpha(k) * Ec; % charlie's Energy in fCB band
sigma_CC = sqrt(lambda * Ec_fCB); 

Charlie_err=0; % Charlie error set to 0
Bob_err=0; % JMAP_error set to 0
total_bits=0; % updates total no. of bits for no. of iterations
p01=0; %probability of charlie decoding alice bit 0 as 1
p10=0; %probability of charlie decoding alice bit 1 as 0

for itr = 1:iterations

% Channel Modelling
hAB = (randn(1, symbols) + 1j * randn(1, symbols)) / sqrt(2);
hCB = (randn(1, symbols) + 1j * randn(1, symbols)) / sqrt(2);
hAC = (sigma_AC / sqrt(2)) * (randn(Nc, symbols) + 1j * randn(Nc, symbols));
hCC = (sigma_CC / sqrt(2)) * (randn(Nc, symbols) + 1j * randn(Nc, symbols));
% Noise Modelling
nB1 = sqrt(N0/2) * (randn(1, symbols)+ 1j * randn(1, symbols));
nB2 = sqrt(N0/2) * (randn(1, symbols) + 1j * randn(1, symbols));
nC1 = sqrt(N0/2) * (randn(Nc, symbols) + 1j * randn(Nc, symbols));

% Initializing Arrays
yB1=zeros(1,symbols);
yB2=zeros(1,symbols);
Bob_x=zeros(1,symbols);
Bob_z=zeros(1,symbols);


% Generating Alice bits
  x_bit = randi([0, 1],1,symbols); 

% Generating Charlie bits
  z_bit = randi([0, M-1],1,symbols);

% TIMESLOT-1
    
    % ALice OOK modulation 
    x_alice = sqrt((1 - alpha(k))*Ec) * x_bit;
    
    % Charlie MQAM modulation 
    z_mod = qammod(z_bit, M, 'UnitAveragePower', true);
    z_charlie = sqrt(alpha(k)*Ec)*z_mod ;

for i= 1 : symbols

    % Rxd at Charlie
     yC1 = x_alice(i) * hAC(:, i) + hCC(:, i) + nC1(:, i);

    % Decoding Error at Charlie    
    x_cap = sum(abs(yC1).^2) > Threshold(k); % Energy Detection at Charlie
    
     % Computing error probabilities
    
    % if x_cap ~= x_bit(i)   % Charlie Errors 
    % 
    % Charlie_err = Charlie_err+1;
    % end
    if x_bit(i)==0 && x_cap==1 
           p01=p01+1;         % Charlie P01 
    elseif x_bit(i)==1 && x_cap==0
           p10=p10+1;        % Charlie P10 
    end

% TIMESLOT - 2

% Modifying Charlie symbols based on decoded Alice symbols
    if x_cap==0
        z_charlie_2 = sqrt((2 - alpha(k)) * Ec) * z_mod(i) * exp(1j * pi/M); % Charlie symbols when alice_cap=0
    else 
        z_charlie_2 = z_mod(i); % Charlie symbols when alice_cap=1
    end  
    
 % Received symbol by Bob
  yB1(i) = hAB(i) * x_alice(i) + nB1(i);
  yB2(i) = hCB(i) * z_charlie_2 + nB2(i);

% JMAP decoder
f11 = (1/(pi * N1b)) * exp(-(abs(yB1(i))^2) / (N1b)); % Conditional PDF of yB1 given x=1
f10 = (1/(pi * N0b)) * exp(-(abs(yB1(i))^2) / (N0b)); % Conditional PDF of yB1 given x=0

f20 = zeros(1, M);  % Conditional PDF of yB2 given x_cap=0
f21 = zeros(1, M);   % Conditional PDF of yB2 given x_cap=1


%Computation of f20, f21 for all M charlie symbols

for j=1:M

z_cap = qammod(j-1, M, 'UnitAveragePower', true);

f21(j) = (1/(pi * N0b)) * exp(-(abs(yB2(i) - sqrt(Ec) * hCB(i) * z_cap)^2) / (N0b));
f20(j) = (1/(pi * N0b)) * exp(-(abs(yB2(i) -sqrt((2 - alpha(k)) * Ec) *hCB(i)* exp(1j * pi/M)*z_cap)^2) / (N0b));
end


% Compute Probabilities 
total_bits = total_bits + symbols;
prob_01=p01/total_bits;
prob_10=p10/total_bits;
prob_00=1-prob_01;
prob_11=1-prob_10;

%JMAP Decoder

fJD0 = f10.*(prob_00 * f20 + prob_01 * f21);
fJD1 = f11.*(prob_10* f20 + prob_11 * f21);
 Charlie_symbol = [fJD0; fJD1]; 
 [~, max_value]= max(Charlie_symbol(:));% Decoding Chalrie and Alice symbol jointly
 [row, col] = ind2sub(size(Charlie_symbol),max_value);
  x_cap_Bob = row -1; % Bob decoded Alice symbol
  z_cap_Bob = col-1;  % Bob decoded Charlie symbol
 
  if x_cap_Bob ~= x_bit(i) || z_cap_Bob ~= z_bit(i)  %JMAP Errors
       Bob_err=Bob_err+1;
  end

end
end
 
% Charlie_p01(Nc,k)=compute_BER(total_bits,p01);%Computing p01 rate
% Charlie_p01(Nc,k)=compute_BER(total_bits,p10);
% Charlie_BER(Nc,k)=compute_BER(total_bits,Charlie_err);
Bob_BER(Nc,k)=compute_BER(total_bits,Bob_err);
   
end
end

figure;
hold on;
C={'mo-','r*-','b+-','yx-'};
for p=1:N

semilogy(alpha, Bob_BER(p,:), C{p}, 'LineWidth', 1, 'MarkerSize', 2);

end
hold off;
xlabel('\alpha'); ylabel('Prob of error ');
title(sprintf('Bob BER versus alpha for SNR = %d dB',snr_dB));
legend('Nc_1', 'Nc_2','Nc_3', 'Nc_4');
grid on;

filename="RESULTS_SNR_35_Nc_1234_itr_1000.mat";
save(filename)
