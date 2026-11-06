function [E,E_P,E_M,E_S] = EEARTH2(TH1,TH2,w_p,w_m,w_s,E_mP,E_mM,E_mS,K_EP,K_EM,K_ES)

% w_p  = the weight factor of the phase 
% w_m  = the weight factor of the magnitude 
% w_s  = the weight factor of the slope  
% E_mP = the maximum allowable percentage of time shift
% E_mM = the maximum allowable magnitude error 
% E_mS = the maximum allowable slope error 
% K_EP = the order of the regression time shift
% K_EM = the order ot the regression magnitude
% K_ES = the order ot the regression slope

% E_P  = EEARTH phase score
% E_M  = EEARTH magnitude score
% E_S  = EEARTH slope score

if nargin == 2
  w_p  = 0.2;
  w_m  = 0.4;
  w_s  = 0.4;
  E_mP = 1;
  E_mM = 1;
  E_mS = 1;
  K_EP = 1;
  K_EM = 1;
  K_ES = 1;
elseif nargin == 5
  E_mP = 1;
  E_mM = 1;
  E_mS = 1;
  K_EP = 1;
  K_EM = 1;
  K_ES = 1;
end

% -------------------------------------------- normalize weighting
w_sum = w_p + w_m + w_s;
w_p = w_p/w_sum;
w_m = w_m/w_sum;
w_s = w_s/w_sum;

% -------------------------------------------- dimension tests
n_t1 = size(TH1,1);
n_t2 = size(TH2,1);
if n_t1 ~= n_t2
  error('** Time histories must have same length !')
end
n_t = n_t1;

n_ddl1 = size(TH1,2);
n_ddl2 = size(TH2,2);
if n_ddl1 ~= n_ddl2
  error('** Time histories must have same number of DOFs !')
end
n_ddl = n_ddl1;

% -------------------------------------------- EEARTH Metric

% EL = number of time shifting steps
% EM = error of magnitude
% ES = error of slope

EL = zeros(1,n_ddl);
EM = zeros(1,n_ddl);
ES = zeros(1,n_ddl);
E_P = zeros(1,n_ddl);
E_M = zeros(1,n_ddl);
E_S = zeros(1,n_ddl);
E = zeros(1,n_ddl);

for i=1:n_ddl

% ------------------------------------------------------------ ith col
  TH1_i = TH1(:,i);
  TH2_i = TH2(:,i);
% ------------------------------------------------------------ Shift signal
  [acor,lags]=xcorr(TH1_i,TH2_i);
  [~,I] = max(abs(acor));
  lag = lags(I);
  if lag < 0
    lag = -lag;
    TH2_i(1:lag) = [];
    TH1_i(end-lag+1:end) = [];
  elseif lag > 0
    TH2_i(end-lag+1:end) = [];
    TH1_i(1:lag) = [];
  end

% ------------------------------------------------------------Error of lags
  EL(i)=abs(lag);
  
% ----------------------------------------------------------- Calculate derivative
  gradth1 = gradient(TH1_i);
  gradth2 = gradient(TH2_i);

% ----------------------------------------------------------- Calculate DTW
  [~,B,C]=dtw(TH1_i,TH2_i);
  [~,G,H]=dtw(gradth1,gradth2);
 
% -------------------------------------------------------Error of magnitude
  if norm(TH2(C),1) == 0
    TH2_i(C) = eps;
    EM(i)=norm((TH1_i(B)-TH2_i(C)),1)/norm(TH2_i(C),1);
  else
    EM(i)=norm((TH1_i(B)-TH2_i(C)),1)/norm(TH2_i(C),1);
  end
  
% -----------------------------------------------------------Error of slope
  if norm(gradth2(H),1) == 0
    gradth2(H) = eps;
     ES(i)=norm((gradth1(G)-gradth2(H)),1)/norm(gradth2(H),1);
  else
     ES(i)= norm((gradth1(G)-gradth2(H)),1)/norm(gradth2(H),1);
  end
    
  if EL(i) == 0
      E_P(i) = 1;
  elseif EL(i) > E_mP*n_t
      E_P(i) = 0;
  else
      E_P(i) = ((E_mP*n_t-EL(i))/(E_mP*n_t))^(K_EP);
  end 

  if EM(i) == 0
      E_M(i) = 1;
  elseif EM(i) > E_mM 
      E_M(i) = 0;
  else
      E_M(i) = ((E_mM-EM(i))/(E_mM))^(K_EM);
  end

  if ES(i) == 0
      E_S(i) = 1;
  elseif ES(i) > E_mS 
      E_S(i) = 0;
  else
      E_S(i) = ((E_mS-ES(i))/(E_mS))^(K_ES);
  end

  % ------------------------------------------ indicator / DOF
  E(i) = w_p*E_P(i) + w_m*E_M(i) + w_s*E_S(i);
end
% -------------------------------------------- mean EEARTH
% E = mean(E);


