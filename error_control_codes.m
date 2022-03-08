% clc;
% clear variables;
% close all;
% 
% N = 1e7;
% n = 7;
% k = 4;
% 
% Nb = N*k;
% 
% P = [1 1 0; 0 1 1; 1 1 1; 1 0 1];
% G = [P eye(k)];
% H = [eye(n-k) P'];
% E = [0 0 0 0 0 0 0;
%      1 0 0 0 0 0 0;
%      0 1 0 0 0 0 0;
%      0 0 1 0 0 0 0;
%      0 0 0 1 0 0 0;
%      0 0 0 0 1 0 0;
%      0 0 0 0 0 1 0;
%      0 0 0 0 0 0 1];
%  
% syn_array = E*H';
% %ms = [0 0 0 0; 0 0 0 1; 0 0 1 0; 0 0 1 1;
% %        0 1 0 0; 0 1 0 1; 0 1 1 0; 0 1 1 1;
% %        1 0 0 0; 1 0 0 1; 1 0 1 0; 1 0 1 1;
% %        1 1 0 0; 1 1 0 1; 1 1 1 0; 1 1 1 1;]
%  
% %cs = mod(ms*G,2);
%    
% %csig = cs;
% %csig(csig==0) =-1;
% syn_val = dec2bin(syn_array');
% syn_val = reshape(syn_val,n-k,2^(n-k));
% syn_dec = bin2dec(syn_val');
% 
% b = randi([0 1], Nb, 1);
% m = reshape(b,k,N)';
% x = mod(m*G,2)';
% s = reshape(x,N*n,1);
% s(s==0)=-1;
% 
% EbNodB = 0:1:15;
% BER0 = 0:1:15;
% BER_Synd = 0:1:15;
% for ii=1:length(EbNodB)
% %     EbNodB(ii)
%     SNRdB = EbNodB(ii)+3.01;
%     rn = awgn(s, SNRdB, 'measured');
%     xcat = rn;
%     xcat(xcat>=0)=1;
%     xcat(xcat<0)=0;
%     xcat = reshape(xcat, n, N);
%     BER0(ii)=sum(sum(xcat~=x))/(N*n);
%     clear xcat;
%     
%     SNRdB = EbNodB(ii)+3.01+10*log10(4/7);
%     rn = awgn(s, SNRdB, 'measured');
%     r = reshape(rn,n,N);
%     y = r;
%     y(y>=0)=1;
%     y(y<0)=0;
%     mH=[]; mE=[];
%     
%     syn=mod(y'*transpose(H),2);
%     syn_val = dec2bin(syn');
%     syn_val = reshape(syn_val, n-k,N);
%     syn_val = bin2dec(syn_val');
%     
%     for jj=1:N
%         e_index(jj)=find(syn_val(jj)==syn_dec);
%     end
%     
%    e = E(e_index,:);
%    scat = mod(y+e',2);
%    bcat = scat(n-k+1:end,:);
%    bcat = reshape(bcat,Nb,1);
%    BER_Synd(ii) = sum(bcat~=b)/Nb;
%    
% end
% 
% %theoretical BER
% 
% EbNo = 10.^(EbNodB/10);
% BER = 1/2*erfc(sqrt(EbNo));
% 
% semilogy(EbNodB,BER,'r*');
% hold on;
% grid on;
% semilogy(EbNodB,BER0,'b-');
% semilogy(EbNodB,BER_Synd,'m-');
% 
% xlabel(' Eb/No (dB)'); ylabel('Bit Error Probablility (P_e)');
% title('BER Performance Statistic');
% axis([0 15 10e-8 1]);
% legend('No Coder(Theoretical)','No Coder(Practical)','Coded');

clc;
clear variables;
close all;

n=7;k=4;

N = 2e5;
Nb = k*N;

P = [1 1 0; 0 1 1; 1 1 1; 1 0 1];
G = [P eye(k)];
H = [eye(n-k) P'];
E = [zeros(1,n);eye(n)];
syn_collection = mod(E*H',2);

bk = randi([0 1],1,Nb);

m = reshape(bk,k,[])';

%Without coding
m_lineCoded = m;
m_lineCoded(m_lineCoded==0)=-1;

bk_lineCoded_woc = reshape(m_lineCoded',1,[]);


%With coding
c = mod(m*G,2);
c_lineCoded = c;
c_lineCoded(c_lineCoded==0) = -1;

bk_lineCoded_wc = reshape(c_lineCoded',1,[]);


SNRdB = -10:1:5;
% SNRdB = 5;
Eb_NodB = SNRdB-10*log10(2);
BER_woc = zeros(1,length(SNRdB));
BER_wc = zeros(1,length(SNRdB));

for i=1:length(SNRdB)
    
    %%%%%%%%%%%%%%%%%Without Coding%%%%%%%%%%%%%%%%%%%
    bk_noisy = awgn(bk_lineCoded_woc,SNRdB(i),'measured');
     
    bk_estimated = bk_noisy;
    bk_estimated(bk_estimated>0) = 1;
    bk_estimated(bk_estimated<0) = 0;
    
    BER_woc(i)=length(find(bk_estimated~=bk))/Nb;
    
    %%%%%%%%%%%%%%%%%With Coding%%%%%%%%%%%%%%%%%%%
    bk_noisy = awgn(bk_lineCoded_wc,SNRdB(i),'measured');
    
    bk_estimated = bk_noisy;
    bk_estimated(bk_estimated>0) = 1;
    bk_estimated(bk_estimated<0) = 0;
    
    c_estimated = reshape(bk_estimated,7,[])';
    syndrome_estimated = mod(c_estimated*H',2);
    
    [a,b]=ismember(syndrome_estimated,syn_collection,'rows');
    E_estimate = E(b,:); 
    
    c_decoded = mod(E_estimate + c_estimated,2);
    m_decoded = c_decoded(:,4:7);
    bk_decoded = reshape(m_decoded',1,[]);
    
    BER_wc(i)=length(find(bk_decoded~=bk))/Nb;
end
semilogy(Eb_NodB,BER_woc,'r:*',Eb_NodB,BER_wc,'b:*');
legend('Without ECC','With ECC');
grid on;
xlabel('Eb/No (dB)');ylabel('BER');