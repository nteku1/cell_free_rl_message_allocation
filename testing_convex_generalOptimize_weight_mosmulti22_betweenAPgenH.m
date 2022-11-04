clear all;


%% Define simulation setup
%processing capabilities
B = 20*10^6; %bandwidth
f =  1*10^6;
C = 20;
%%%D = 10*10^6;
S_u= 1*10^6;

%Number of Monte Carlo setups
nbrOfSetups = 4;%3 30;
aaaaaa = 0;
outputs = zeros(1,nbrOfSetups);
indices = zeros(1,nbrOfSetups);
%Number of channel realizations per setup
nbrOfRealizations = 1; %1000;
%Number of APs in the cell-free network
L = 10;%10; %original: 10

%Number of UEs
%K = 2;
K = 3;

%Number of antennas per AP
N = 1;

%Length of the coherence block
tau_c = 200;

%Number of pilots per coherence block
tau_p = 20;

%Uplink transmit power per UE (mW)
p = 100;

%power ratio
ada_1 = 1;
beta1 = 10;
beta2 = 1;%0.5;
% beta3 = 10^3;
% beta4 = 10^7;
% beta5 = 10^7;
% beta3 = 10^3;
% beta4 = 10^7;
% fileID3 = fopen('VER1_10AP_3UE_NO_SHAD_Hmat_scenario3_10APs_Multi_2_users_FUNFINALACTUREG_complex_50000_7_28_ExtendingTesting_ver62.txt','w');
% fileID4 = fopen('VER1_10AP_3UE_NO_SHAD_Hmat_scenario3_10APs_Multi_2_users_FUNFINALACTUREG_complex_part2_7_28_ExtendingTesting_ver62.txt','w');

% fileID3 = fopen('VER3_RRRREDO_SHADDD_Hmat_scenario3_10APs_Multi_3_users_FUNFINALACTUREG_complex_50000_7_28_ExtendingTesting_ver62.txt','w');
% fileID4 = fopen('VER3_RRRREDO_SHADDD_Hmat_scenario3_10APs_Multi_3_users_FUNFINALACTUREG_complex_part2_7_28_ExtendingTesting_ver62.txt','w');
% fileID2 = fopen('alphasssnew.txt','w');
%log2(1+((p*ada_1*alphasss^2*abs(H_AP).^2)/(1+p*ada_1*((1-alphasss)^2)*abs(H_AP)^2)))
for iiii = 1:51000%210000
    if iiii == 1 %initialize both user and AP positions
        [gainOverNoisedB,R,pilotIndexCF,pilotIndexSC,APpositions,UEpositions] = generateSetup_threeslope_rev(L,K,N,tau_p,1,p);
    else
         [gainOverNoisedB,R,pilotIndexCF,pilotIndexSC] = generateSetup_threeslope_rev_justuserpos_change22(L,K,N,tau_p,1,p,APpositions,UEpositions); 
    end
        betaVal = db2pow(gainOverNoisedB);
   [Hhat_AP,H_AP,B_AP] = functionChannelEstimates(R,nbrOfRealizations,L,K,N,tau_p,pilotIndexCF,p);
%     for n = 1:nbrOfSetups
%         SINR = 1:0.11*n:n+1;
%         indices(n) = max(SINR)-min(SINR);
%         D = 1 *10^6;
% 
% %        [gainOverNoisedB,R,pilotIndexCF,pilotIndexSC] = generateSetup_threeslope(L,K,N,tau_p,1,p);
% %         betaVal = db2pow(gainOverNoisedB);
% 
% 
%         %Full transmit power case
% 
%         %Generate channel realizations, channel estimates, and estimation
%         %error correlation matrices for all UEs to the APs 
% %         [Hhat_AP,H_AP,B_AP] = functionChannelEstimates(R,nbrOfRealizations,L,K,N,tau_p,pilotIndexCF,p);
%         %old
%         %fun = @(alphasss) max(((alphasss*D*C)/f) + ((alphasss*S_u)/(B*log2(1+((p*ada_1*alphasss.^2*abs(H_AP).^2)/(1+p*ada_1*((1-alphasss).^2)*abs(H_AP).^2))))));
%         %initial condition
%        
% %         fun = @(alphasss) ((p*ada_1*alphasss.^2*abs(H_AP).^2))/(1+p*ada_1*((1-alphasss).^2)*abs(H_AP).^2);
% %         fun = @(alphasss) sum(((log2(1+((p*ada_1*alphasss.^2*abs(H_AP).^2)/(1+p*ada_1*((1-alphasss).^2)*abs(H_AP).^2))))));
%        %max
%         
%          %added latency of sending H_AP - may need to revise rate for this?
%          
%          %Hhat_AP + sending fragement from AP to central AP + sending
%          %alphas!!!
%          %(beta3*(sum(alphasss.*D)/sum((B*log2(1+((p*ada_1*alphasss.^2.*abs(H_AP).^2.')./(1+p*ada_1*((1-alphasss).^2).*abs(H_AP).^2.'))))))) 
% %            fun = @(alphasss) max((beta1*((alphasss.*D*C)/f)) + (beta4*(((L-1)*(length(H_AP)+length(alphasss)))/sum((B*log2(1+((p*ada_1*alphasss.^2.*abs(H_AP).^2.')./(1+p*ada_1*((1-alphasss).^2).*abs(H_AP).^2.')))))))  +   (beta5*(((L-1)*(length(alphasss)))/sum((B*log2(1+((p*ada_1*alphasss.^2.*abs(H_AP).^2.')./(1+p*ada_1*((1-alphasss).^2).*abs(H_AP).^2.'))))))) + (beta2*((alphasss*S_u)./(B*log2(1+((p*ada_1*alphasss.^2.*abs(H_AP).^2.')./(1+p*ada_1*((1-alphasss).^2).*abs(H_AP).^2.')))))));    
%           fun = @(alphasss) max((beta1*((alphasss.*D*C)/f)) + (beta2*((alphasss*S_u)./(B*log2(1+((p*ada_1*alphasss.^2.*abs(H_AP).^2.')./(1+p*ada_1*((1-alphasss).^2).*abs(H_AP).^2.')))))));    
%                   
% 
% 
% %          fun = @(alphasss) max((beta1*((alphasss.*D*C)/f)) + (beta3*((alphasss.*(L-1)*D)/sum((B*log2(1+((p*ada_1*alphasss.^2.*abs(Hhat_AP).^2.')./(1+p*ada_1*((1-alphasss).^2).*abs(Hhat_AP).^2.')))))))  + (beta4*(((L-1)*(length(Hhat_AP)+length(alphasss)))/sum((B*log2(1+((p*ada_1*alphasss.^2.*abs(Hhat_AP).^2.')./(1+p*ada_1*((1-alphasss).^2).*abs(Hhat_AP).^2.'))))))) + (beta2*((alphasss*S_u)./(B*log2(1+((p*ada_1*alphasss.^2.*abs(Hhat_AP).^2.')./(1+p*ada_1*((1-alphasss).^2).*abs(Hhat_AP).^2.')))))));
%          
%          
% 
% 
% 
% 
% 
% 
% 
% 
%          %Replace H_AP with Hhat_AP
% %                            fun = @(alphasss) max((beta1*((alphasss.*D*C)/f)) + (((L-1)*length(Hhat_AP))/sum((B*log2(1+((p*ada_1*alphasss.^2.*abs(Hhat_AP).^2.')./(1+p*ada_1*((1-alphasss).^2).*abs(Hhat_AP).^2.')))))) + (beta2*((alphasss*S_u)./(B*log2(1+((p*ada_1*alphasss.^2.*abs(Hhat_AP).^2.')./(1+p*ada_1*((1-alphasss).^2).*abs(Hhat_AP).^2.')))))));
%          
%          %Replace Length(H_AP) with product of dimensions of H_AP or * 2?
% %                   fun = @(alphasss) max((beta1*((alphasss.*D*C)/f)) + (((L-1)*length(H_AP))/sum((B*log2(1+((p*ada_1*alphasss.^2.*abs(H_AP).^2.')./(1+p*ada_1*((1-alphasss).^2).*abs(H_AP).^2.')))))) + (beta2*((alphasss*S_u)./(B*log2(1+((p*ada_1*alphasss.^2.*abs(H_AP).^2.')./(1+p*ada_1*((1-alphasss).^2).*abs(H_AP).^2.')))))));
%          
%           %added latency of AP to AP sending decoded fragments
% %          fun = @(alphasss) max((beta1*((alphasss.*D*C)/f)) + ((alphasss.*(L-1)*D)/sum((B*log2(1+((p*ada_1*alphasss.^2.*abs(H_AP).^2.')./(1+p*ada_1*((1-alphasss).^2).*abs(H_AP).^2.')))))) + (beta2*((alphasss*S_u)./(B*log2(1+((p*ada_1*alphasss.^2.*abs(H_AP).^2.')./(1+p*ada_1*((1-alphasss).^2).*abs(H_AP).^2.')))))));
%          
%          
% % % % % %         fun = @(alphasss) max(((alphasss.*D*C)/f) + ((alphasss*S_u)./(B*log2(1+((p*ada_1*alphasss.^2.*abs(H_AP).^2.')./(1+p*ada_1*((1-alphasss).^2).*abs(H_AP).^2.'))))));
% % % % % %         fun = @(alphasss) max(((alphasss.*D*C)/f) + ((alphasss*S_u)./(B*log2(1+SINR)))); 
%         if n == 1
%         alphasss_0= rand(1,L);
%         alphasss_0= alphasss_0/sum(alphasss_0);
%         outputs(n) = outputs(n) + fun(alphasss_0);
% 
%         else
%         alphasss_0 = alphasss;
%         end
% 
%         lb = zeros(1,L);%0;
%         ub = ones(1,L);%;1;
%         Aeq = ones(1,L);
%         beq = 1;
%         A = [];
%         b = [];
% 
%         if n > 1
%             outputs(n) = outputs(n) + fun(alphasss);
%         end
% % % %        alphasss = fmincon(fun,alphasss_0,A,b,Aeq,beq,lb,ub);
% % % %        %'MaxFunEvals',4500,'TolCon', 1e-10'  6500
% 
% 
%       options = optimoptions(@fmincon,'MaxFunEvals',4500,'StepTolerance',1e-15,'TolCon',1e-10);%'MaxIter',2000,'TolCon',1e-10); %2000 TolCon 1e-10 TolX
%        alphasss = fmincon(fun,alphasss_0,A,b,Aeq,beq,lb,ub,[],options);
%  
% 
%         d= 1;
%         
%        
% 
% 
%     end
%     
     %%%%New Stuff put Hhat or no

     
%      SINRssss = ((p*ada_1*alphasss.^2.*abs(H_AP).^2.')./(1+p*ada_1*((1-alphasss).^2).*abs(H_AP).^2.'));
%         [IIIIIII,BBBBBBB1] = sort(SINRssss,'ascend');
%         [IIIIIII2,BBBBBBB2] = sort(alphasss,'ascend');
%         
%         if any(BBBBBBB1 ~= BBBBBBB2)
%             aaaaaa = aaaaaa + 1;
%         end
%         fprintf(fileID,'%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f\n',SINRssss(1),...
%             SINRssss(2),SINRssss(3),SINRssss(4),SINRssss(5),SINRssss(6),SINRssss(7),...
%             SINRssss(8),SINRssss(9),SINRssss(10));
%         
%          fprintf(fileID2,'%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f\n',alphasss(1),...
%             alphasss(2),alphasss(3),alphasss(4),alphasss(5),alphasss(6),alphasss(7),...
%             alphasss(8),alphasss(9),alphasss(10));
%     
    
if iiii <= 50000% 1000
% fprintf(fileID3,'%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f\n',...
%          real(H_AP(1)), imag(H_AP(1)),real(H_AP(2)),imag(H_AP(2)),real(H_AP(3)), imag(H_AP(3)),real(H_AP(4)), imag(H_AP(4)),real(H_AP(5)), imag(H_AP(5)));%,real(H_AP(6)), imag(H_AP(6)),real(H_AP(7)), imag(H_AP(7)));    

% % % % fprintf(fileID3,'%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f\n',...
% % % %           real(H_AP(1)), imag(H_AP(1)),real(H_AP(2)),imag(H_AP(2)),real(H_AP(3)), imag(H_AP(3)),...
% % % %           real(H_AP(4)), imag(H_AP(4)),real(H_AP(5)),imag(H_AP(5)),real(H_AP(6)), imag(H_AP(6)),...
% % % %           real(H_AP(7)), imag(H_AP(7)),real(H_AP(8)),imag(H_AP(8)),real(H_AP(9)), imag(H_AP(9)),...
% % % %           real(H_AP(10)), imag(H_AP(10)));

for jlock = 1:K
fprintf(fileID3,'%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f\n',...
          real(H_AP(1,jlock)), imag(H_AP(1,jlock)),real(H_AP(2,jlock)),imag(H_AP(2,jlock)),real(H_AP(3,jlock)), imag(H_AP(3,jlock)),...
          real(H_AP(4,jlock)), imag(H_AP(4,jlock)),real(H_AP(5,jlock)),imag(H_AP(5,jlock)),real(H_AP(6,jlock)), imag(H_AP(6,jlock)),...
          real(H_AP(7,jlock)), imag(H_AP(7,jlock)),real(H_AP(8,jlock)),imag(H_AP(8,jlock)),real(H_AP(9,jlock)), imag(H_AP(9,jlock)),...
          real(H_AP(10,jlock)), imag(H_AP(10,jlock)));

 % 15 APs
%  fprintf(fileID3,'%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f\n',...
%   real(H_AP(1,jlock)), imag(H_AP(1,jlock)),real(H_AP(2,jlock)),imag(H_AP(2,jlock)),real(H_AP(3,jlock)), imag(H_AP(3,jlock)),...
%   real(H_AP(4,jlock)), imag(H_AP(4,jlock)),real(H_AP(5,jlock)),imag(H_AP(5,jlock)),real(H_AP(6,jlock)), imag(H_AP(6,jlock)),...
%   real(H_AP(7,jlock)), imag(H_AP(7,jlock)),real(H_AP(8,jlock)),imag(H_AP(8,jlock)),real(H_AP(9,jlock)), imag(H_AP(9,jlock)),...
%   real(H_AP(10,jlock)), imag(H_AP(10,jlock)),real(H_AP(11,jlock)), imag(H_AP(11,jlock)),real(H_AP(12,jlock)), imag(H_AP(12,jlock)),...
%   real(H_AP(13,jlock)), imag(H_AP(13,jlock)),real(H_AP(14,jlock)), imag(H_AP(14,jlock)),real(H_AP(15,jlock)), imag(H_AP(15,jlock)));  
%      
%       
      
      
% fprintf(fileID3,'%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f\n',...
%   real(H_AP(1,jlock)), imag(H_AP(1,jlock)),real(H_AP(2,jlock)),imag(H_AP(2,jlock)),real(H_AP(3,jlock)), imag(H_AP(3,jlock)),...
%   real(H_AP(4,jlock)), imag(H_AP(4,jlock)),real(H_AP(5,jlock)),imag(H_AP(5,jlock)),real(H_AP(6,jlock)), imag(H_AP(6,jlock)),...
%   real(H_AP(7,jlock)), imag(H_AP(7,jlock)),real(H_AP(8,jlock)),imag(H_AP(8,jlock)),real(H_AP(9,jlock)), imag(H_AP(9,jlock)),...
%   real(H_AP(10,jlock)), imag(H_AP(10,jlock)),real(H_AP(11,jlock)), imag(H_AP(11,jlock)),real(H_AP(12,jlock)), imag(H_AP(12,jlock)),...
%   real(H_AP(13,jlock)), imag(H_AP(13,jlock)),real(H_AP(14,jlock)), imag(H_AP(14,jlock)),real(H_AP(15,jlock)), imag(H_AP(15,jlock)),...
%   real(H_AP(16,jlock)), imag(H_AP(16,jlock)),real(H_AP(17,jlock)), imag(H_AP(17,jlock)),real(H_AP(18,jlock)), imag(H_AP(18,jlock)),...
%   real(H_AP(19,jlock)), imag(H_AP(19,jlock)),real(H_AP(20,jlock)), imag(H_AP(20,jlock)));  


% fprintf(fileID3,'%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f\n',...
%   real(H_AP(1,jlock)), imag(H_AP(1,jlock)),real(H_AP(2,jlock)),imag(H_AP(2,jlock)),real(H_AP(3,jlock)), imag(H_AP(3,jlock)),...
%   real(H_AP(4,jlock)), imag(H_AP(4,jlock)),real(H_AP(5,jlock)),imag(H_AP(5,jlock)),real(H_AP(6,jlock)), imag(H_AP(6,jlock)),...
%   real(H_AP(7,jlock)), imag(H_AP(7,jlock)),real(H_AP(8,jlock)),imag(H_AP(8,jlock)),real(H_AP(9,jlock)), imag(H_AP(9,jlock)),...
%   real(H_AP(10,jlock)), imag(H_AP(10,jlock)),real(H_AP(11,jlock)), imag(H_AP(11,jlock)),real(H_AP(12,jlock)), imag(H_AP(12,jlock)),...
%   real(H_AP(13,jlock)), imag(H_AP(13,jlock)),real(H_AP(14,jlock)), imag(H_AP(14,jlock)),real(H_AP(15,jlock)), imag(H_AP(15,jlock)),...
%   real(H_AP(16,jlock)), imag(H_AP(16,jlock)),real(H_AP(17,jlock)), imag(H_AP(17,jlock)),real(H_AP(18,jlock)), imag(H_AP(18,jlock)),...
%   real(H_AP(19,jlock)), imag(H_AP(19,jlock)),real(H_AP(20,jlock)), imag(H_AP(20,jlock)),real(H_AP(21,jlock)), imag(H_AP(21,jlock)),... 
%   real(H_AP(22,jlock)), imag(H_AP(22,jlock)),real(H_AP(23,jlock)), imag(H_AP(23,jlock)),real(H_AP(24,jlock)), imag(H_AP(24,jlock)),...
%   real(H_AP(25,jlock)), imag(H_AP(25,jlock)));  

%5 APs
% fprintf(fileID3,'%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f\n',...
%   real(H_AP(1,jlock)), imag(H_AP(1,jlock)),real(H_AP(2,jlock)),imag(H_AP(2,jlock)),real(H_AP(3,jlock)), imag(H_AP(3,jlock)),...
%   real(H_AP(4,jlock)), imag(H_AP(4,jlock)),real(H_AP(5,jlock)),imag(H_AP(5,jlock)));  
end
      
      
else
% fprintf(fileID4,'%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f\n',...
%          real(H_AP(1)), imag(H_AP(1)),real(H_AP(2)),imag(H_AP(2)),real(H_AP(3)), imag(H_AP(3)),real(H_AP(4)), imag(H_AP(4)),real(H_AP(5)), imag(H_AP(5)));%,real(H_AP(6)), imag(H_AP(6)),real(H_AP(7)), imag(H_AP(7)));  

% % % % % % % % fprintf(fileID4,'%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f\n',...
% % % % % % % %           real(H_AP(1)), imag(H_AP(1)),real(H_AP(2)),imag(H_AP(2)),real(H_AP(3)), imag(H_AP(3)),...
% % % % % % % %           real(H_AP(4)), imag(H_AP(4)),real(H_AP(5)),imag(H_AP(5)),real(H_AP(6)), imag(H_AP(6)),...
% % % % % % % %           real(H_AP(7)), imag(H_AP(7)),real(H_AP(8)),imag(H_AP(8)),real(H_AP(9)), imag(H_AP(9)),...
% % % % % % % %           real(H_AP(10)), imag(H_AP(10)));


for block = 1:K
fprintf(fileID4,'%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f\n',...
          real(H_AP(1,block)), imag(H_AP(1,block)),real(H_AP(2,block)),imag(H_AP(2,block)),real(H_AP(3,block)), imag(H_AP(3,block)),...
          real(H_AP(4,block)), imag(H_AP(4,block)),real(H_AP(5,block)),imag(H_AP(5,block)),real(H_AP(6,block)), imag(H_AP(6,block)),...
          real(H_AP(7,block)), imag(H_AP(7,block)),real(H_AP(8,block)),imag(H_AP(8,block)),real(H_AP(9,block)), imag(H_AP(9,block)),...
          real(H_AP(10,block)), imag(H_AP(10,block)));
   

% fprintf(fileID4,'%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f\n',...
%           real(H_AP(1,block)), imag(H_AP(1,block)),real(H_AP(2,block)),imag(H_AP(2,block)),real(H_AP(3,block)), imag(H_AP(3,block)),...
%           real(H_AP(4,block)), imag(H_AP(4,block)),real(H_AP(5,block)),imag(H_AP(5,block)),real(H_AP(6,block)), imag(H_AP(6,block)),...
%           real(H_AP(7,block)), imag(H_AP(7,block)),real(H_AP(8,block)),imag(H_AP(8,block)),real(H_AP(9,block)), imag(H_AP(9,block)),...
%           real(H_AP(10,block)), imag(H_AP(10,block)),real(H_AP(11,jlock)), imag(H_AP(11,jlock)),real(H_AP(12,jlock)), imag(H_AP(12,jlock)),...
%           real(H_AP(13,jlock)), imag(H_AP(13,jlock)),real(H_AP(14,jlock)), imag(H_AP(14,jlock)),real(H_AP(15,jlock)), imag(H_AP(15,jlock)));     
%       
%       
      
%       
% fprintf(fileID4,'%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f\n',...
%           real(H_AP(1,block)), imag(H_AP(1,block)),real(H_AP(2,block)),imag(H_AP(2,block)),real(H_AP(3,block)), imag(H_AP(3,block)),...
%           real(H_AP(4,block)), imag(H_AP(4,block)),real(H_AP(5,block)),imag(H_AP(5,block)),real(H_AP(6,block)), imag(H_AP(6,block)),...
%           real(H_AP(7,block)), imag(H_AP(7,block)),real(H_AP(8,block)),imag(H_AP(8,block)),real(H_AP(9,block)), imag(H_AP(9,block)),...
%           real(H_AP(10,block)), imag(H_AP(10,block)),real(H_AP(11,jlock)), imag(H_AP(11,jlock)),real(H_AP(12,jlock)), imag(H_AP(12,jlock)),...
%           real(H_AP(13,jlock)), imag(H_AP(13,jlock)),real(H_AP(14,jlock)), imag(H_AP(14,jlock)),real(H_AP(15,jlock)), imag(H_AP(15,jlock)),...
%           real(H_AP(16,jlock)), imag(H_AP(16,jlock)),real(H_AP(17,jlock)), imag(H_AP(17,jlock)),real(H_AP(18,jlock)), imag(H_AP(18,jlock)),...
%           real(H_AP(19,jlock)), imag(H_AP(19,jlock)),real(H_AP(20,jlock)), imag(H_AP(20,jlock)));     
% 

% fprintf(fileID4,'%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f\n',...
%           real(H_AP(1,block)), imag(H_AP(1,block)),real(H_AP(2,block)),imag(H_AP(2,block)),real(H_AP(3,block)), imag(H_AP(3,block)),...
%           real(H_AP(4,block)), imag(H_AP(4,block)),real(H_AP(5,block)),imag(H_AP(5,block)),real(H_AP(6,block)), imag(H_AP(6,block)),...
%           real(H_AP(7,block)), imag(H_AP(7,block)),real(H_AP(8,block)),imag(H_AP(8,block)),real(H_AP(9,block)), imag(H_AP(9,block)),...
%           real(H_AP(10,block)), imag(H_AP(10,block)),real(H_AP(11,jlock)), imag(H_AP(11,jlock)),real(H_AP(12,jlock)), imag(H_AP(12,jlock)),...
%           real(H_AP(13,jlock)), imag(H_AP(13,jlock)),real(H_AP(14,jlock)), imag(H_AP(14,jlock)),real(H_AP(15,jlock)), imag(H_AP(15,jlock)),...
%           real(H_AP(16,jlock)), imag(H_AP(16,jlock)),real(H_AP(17,jlock)), imag(H_AP(17,jlock)),real(H_AP(18,jlock)), imag(H_AP(18,jlock)),...
%           real(H_AP(19,jlock)), imag(H_AP(19,jlock)),real(H_AP(20,jlock)), imag(H_AP(20,jlock)),real(H_AP(21,jlock)), imag(H_AP(21,jlock)),... 
%           real(H_AP(22,jlock)), imag(H_AP(22,jlock)),real(H_AP(23,jlock)), imag(H_AP(23,jlock)),real(H_AP(24,jlock)), imag(H_AP(24,jlock)),...
%           real(H_AP(25,jlock)), imag(H_AP(25,jlock)));     



%5 APs
% fprintf(fileID4,'%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f\n',...
%           real(H_AP(1,block)), imag(H_AP(1,block)),real(H_AP(2,block)),imag(H_AP(2,block)),real(H_AP(3,block)), imag(H_AP(3,block)),...
%           real(H_AP(4,block)), imag(H_AP(4,block)),real(H_AP(5,block)),imag(H_AP(5,block)));      
end



end
    
end
fclose(fileID3);
fclose(fileID4);
% outputs = outputs/iiii
% plot(1:nbrOfSetups,outputs,'-bo');
% xlabel('Timesteps');
% ylabel('Max Latency');
% title('Latency (\mu_1 = 100, \mu_2 = 1000 vs. Timesteps - Scenario 1: APs unchanged'); %\beta =-81.2 dB 3 w/ Shadow Fading
% %title('Optimizing \alpha - Constant H'); \mu_3 = 10^3, \mu_4 = \mu_5 = 10^7) 
% fclose(fileID);
% fclose(fileID2);