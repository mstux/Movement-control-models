clear

tic;

Ne = 80; % excitatory neurons
Ni = 20; % inhibitory neurons
PMv = pop(Ne,Ni);
AIP = pop(Ne,Ni);

proj = .4; % how many projection neurons, same for AIP and PMv
aff = .4; % how many neurons receiving afferents, same for AIP and PMv
w_in_PMv1 = zeros(Ne+Ni,1); % input weights
for i=1:Ne+Ni
    if mod(i,2)==0
        w_in_PMv1(i) = normrnd(2.5,.75);
    end
%     if rand<.4
%         w_in_PMv1(i) = normrnd(2.5,.75);
%     end
end
w_in_PMv2 = zeros(Ne+Ni,1); % input weights
for i=1:Ne+Ni
    if mod(i,2)==1
        w_in_PMv2(i) = normrnd(2.5,.75);
    end
%     if rand<.4
%         w_in_PMv2(i) = normrnd(2.5,.75);
%     end
end
w_in_AIP1 = zeros(Ne+Ni,1); % input weights
for i=1:Ne+Ni
    if mod(i,2)==1
        w_in_AIP1(i) =  normrnd(2.5,.75);
    end
%     if rand<.4
%         w_in_AIP1(i) =  normrnd(2.5,.75);
%     end
end
w_in_AIP2 = zeros(Ne+Ni,1); % input weights
for i=1:Ne+Ni
    if mod(i,2)==0
        w_in_AIP2(i) =  normrnd(2.5,.75);
    end
%     if rand<.4
%         w_in_AIP2(i) =  normrnd(2.5,.75);
%     end
end

% sparse connectivity
mask_PMv = ones(Ne+Ni,Ne+Ni);
mask_AIP = ones(Ne+Ni,Ne+Ni);
for i=1:Ne+Ni
    for j=1:Ne+Ni
        if rand<.75
            mask_PMv(i,j) = 0;
        end
        if rand<.75
            mask_AIP(i,j) = 0;
        end
    end
end
PMv.w(mask_PMv==0) = 0;
AIP.w(mask_AIP==0) = 0;

% inter-area connections
w_PMvtoAIP = normrnd(1.5,.4, (Ne+Ni)*aff, (Ne+Ni)*proj );
w_AIPtoPMv = normrnd(1.5,.4, (Ne+Ni)*aff, (Ne+Ni)*proj );

wmax = 10; % weights upper bound

% short-term plasticity
STD_PMvtoAIP = ones( (Ne+Ni)*aff, (Ne+Ni)*proj );
STD_AIPtoPMv = ones( (Ne+Ni)*aff, (Ne+Ni)*proj );
% PMv.w(1:(Ne+Ni)*proj,(Ne+Ni)*proj+1:(Ne+Ni)*proj+(Ne+Ni)*aff) = 1.1*PMv.w(1:(Ne+Ni)*proj,(Ne+Ni)*proj+1:(Ne+Ni)*proj+(Ne+Ni)*aff);
% AIP.w(1:(Ne+Ni)*proj,(Ne+Ni)*proj+1:(Ne+Ni)*proj+(Ne+Ni)*aff) = 1.1*AIP.w(1:(Ne+Ni)*proj,(Ne+Ni)*proj+1:(Ne+Ni)*proj+(Ne+Ni)*aff);

firings_PMv = [];                      % spike timings 
firings_AIP = []; 

% Plasticity model parameters
Ad = 0.00014;                       % Depression Amplitude
Ap = 0.00008;                       % Potentiation Amplitude
tau_p = 7;                             % time constant of the voltage low-pass for the potentiation term
tau_x = 15;                            % time constant for the presynaptic spike low-pass
tau_d = 10;                            % time constant of the voltage low-pass for the depression term
theta_minus = -70;                     
theta_plus = -45;

% Short term plasticity model parameters
U = .2; % synapse depotentiation fraction after a spike
taud = 750;

for t=1:60000          % simulation of 60 s

  t  
    
  I_AIP=[3*randn(Ne,1); 2*randn(Ni,1)]; % noise
  for s=0:9 % 10 plus one stimulations
      if t>10000+3000*s && t<=10400+3000*s % stimulus for 400 ms
          if mod(s,2)==0
              I_AIP = I_AIP+w_in_AIP1;
          else
              I_AIP = I_AIP+w_in_AIP2;
          end
      end
  end
  if t>55000 && t<=55400 % stimulus for 400 ms
      I_AIP = I_AIP+w_in_AIP2;
  end
  if t>50000 && t<=50400 % stimulus for 400 ms
      I_AIP = I_AIP+w_in_AIP1;
  end
  fired_AIP=find(AIP.v>=30);    % indices of spikes
  firings_AIP=[firings_AIP; t+0*fired_AIP,fired_AIP];
  AIP.v(fired_AIP)=AIP.c(fired_AIP);
  AIP.u(fired_AIP)=AIP.u(fired_AIP)+AIP.d(fired_AIP);
  w_STD_AIP = [AIP.w(:,1:Ne).*AIP.STD, AIP.w(:,Ne+1:end)];
  I_AIP=I_AIP+sum(w_STD_AIP(:,fired_AIP),2);
    
  I_PMv=[3*randn(Ne,1); 2*randn(Ni,1)]; % noise
  for s=0:9 % 10 stimulations
      if t>10000+3000*s && t<=10400+3000*s % stimulus for 400 ms
          if mod(s,2)==0
              I_PMv = I_PMv+w_in_PMv1;
          else
              I_PMv = I_PMv+w_in_PMv2;
          end
      end
  end
  fired_PMv=find(PMv.v>=30);    % indices of spikes
  firings_PMv=[firings_PMv; t+0*fired_PMv,fired_PMv];
  PMv.v(fired_PMv)=PMv.c(fired_PMv);
  PMv.u(fired_PMv)=PMv.u(fired_PMv)+PMv.d(fired_PMv);
  w_STD_PMv = [PMv.w(:,1:Ne).*PMv.STD, PMv.w(:,Ne+1:end)];
  I_PMv=I_PMv+sum(w_STD_PMv(:,fired_PMv),2);
  
  STD_w_AIPtoPMv = STD_AIPtoPMv.*w_AIPtoPMv;
  STD_w_PMvtoAIP = STD_PMvtoAIP.*w_PMvtoAIP;
  I_PMv(1+(Ne+Ni)*proj:(Ne+Ni)*proj+(Ne+Ni)*aff)=I_PMv(1+(Ne+Ni)*proj:(Ne+Ni)*proj+(Ne+Ni)*aff)+sum(STD_w_AIPtoPMv( :, fired_AIP( find(fired_AIP<=(Ne+Ni)*proj) )),2); % if for example proj=aff=.2 and Ne+Ni=100 neurons 21:40 of PMv receive from neurons 1:20 of AIP
  I_AIP(1+(Ne+Ni)*proj:(Ne+Ni)*proj+(Ne+Ni)*aff)=I_AIP(1+(Ne+Ni)*proj:(Ne+Ni)*proj+(Ne+Ni)*aff)+sum(STD_w_PMvtoAIP( :, fired_PMv( find(fired_PMv<=(Ne+Ni)*proj) )),2); % if for example proj=aff=.2 and Ne+Ni=100 neurons 21:40 of AIP receive from neurons 1:20 of PMv
  
%   I_PMv(1) = I_PMv(1) + 3.5*randn;
%   I_AIP(21) = I_AIP(21) + 3.5*randn;
  
  % step 0.5 ms for numerical stability PMv
  PMv.v = PMv.v+0.5*(0.04*PMv.v.^2+5*PMv.v+140-PMv.u+I_PMv); 
  PMv.v = PMv.v+0.5*(0.04*PMv.v.^2+5*PMv.v+140-PMv.u+I_PMv);
  reset_PMv=find(PMv.v>=30);
  PMv.v(reset_PMv)=30; % standardize action potentials
  PMv.u=PMv.u+PMv.a.*(PMv.b.*PMv.v-PMv.u);
  
  % step 0.5 ms for numerical stability AIP
  AIP.v = AIP.v+0.5*(0.04*AIP.v.^2+5*AIP.v+140-AIP.u+I_AIP); 
  AIP.v = AIP.v+0.5*(0.04*AIP.v.^2+5*AIP.v+140-AIP.u+I_AIP);
  reset_AIP=find(AIP.v>=30);
  AIP.v(reset_AIP)=30; % standardize action potentials
  AIP.u=AIP.u+AIP.a.*(AIP.b.*AIP.v-AIP.u);

  % plasticity PMv
  fired_e_PMv = fired_PMv(fired_PMv<=Ne); % excitatory neurons that fired
  PMv.v_md = PMv.v_md+(PMv.v-PMv.v_md)/tau_d;
  PMv.v_mp = PMv.v_mp+(PMv.v-PMv.v_mp)/tau_p;
  PMv.v_rekt_md = ((PMv.v_md-theta_minus) > 0).*(PMv.v_md-theta_minus); % threshold the low pass
  PMv.v_rekt_mp = ((PMv.v_mp-theta_minus) > 0).*(PMv.v_mp-theta_minus); % threshold the low pass
  PMv.vy = ((PMv.v-theta_plus) > 0).*(PMv.v-theta_plus); % threshold of voltage
  PMv.x0(:,fired_e_PMv) = 1;
  PMv.x = PMv.x+(PMv.x0-PMv.x)/tau_x;
  PMv.w(:,1:Ne) = PMv.w(:,1:Ne)-Ad*PMv.x0.*PMv.v_rekt_md+Ap*PMv.x.*PMv.vy.*PMv.v_rekt_mp; % weight update
  idx_PMv = find(PMv.w(:,1:Ne)<0); % weight lower bound
  w_temp_PMv = PMv.w(:,1:Ne); 
  w_temp_PMv(idx_PMv) = 0;
  PMv.w(:,1:Ne) = w_temp_PMv;
  PMv.w(PMv.w>PMv.wmax) = PMv.wmax; % weight upper bound

  % plasticity AIP
  fired_e_AIP = fired_AIP(fired_AIP<=Ne); % excitatory neurons that fired
  AIP.v_md = AIP.v_md+(AIP.v-AIP.v_md)/tau_d;
  AIP.v_mp = AIP.v_mp+(AIP.v-AIP.v_mp)/tau_p;
  AIP.v_rekt_md = ((AIP.v_md-theta_minus) > 0).*(AIP.v_md-theta_minus); % threshold the low pass
  AIP.v_rekt_mp = ((AIP.v_mp-theta_minus) > 0).*(AIP.v_mp-theta_minus); % threshold the low pass
  AIP.vy = ((AIP.v-theta_plus) > 0).*(AIP.v-theta_plus); % threshold of voltage
  AIP.x0(:,fired_e_AIP) = 1;
  AIP.x = AIP.x+(AIP.x0-AIP.x)/tau_x;
  AIP.w(:,1:Ne) = AIP.w(:,1:Ne)-Ad*AIP.x0.*AIP.v_rekt_md+Ap*AIP.x.*AIP.vy.*AIP.v_rekt_mp; % weight update
  idx_AIP = find(AIP.w(:,1:Ne)<0); % weight lower bound
  w_temp_AIP = AIP.w(:,1:Ne); 
  w_temp_AIP(idx_AIP) = 0;
  AIP.w(:,1:Ne) = w_temp_AIP;
  AIP.w(AIP.w>AIP.wmax) = AIP.wmax; % weight upper bound
  
  % STP PMv
  PMv.STD = PMv.STD +(1-PMv.STD)/taud -U*PMv.x0;
  PMv.STD(PMv.STD<0) = 0;
  
  % STP AIP
  AIP.STD = AIP.STD +(1-AIP.STD)/taud -U*AIP.x0;
  AIP.STD(AIP.STD<0) = 0;
  
%   % STP AIPtoPMv
%   STD_AIPtoPMv = STD_AIPtoPMv +(1-STD_AIPtoPMv)/taud -U*AIP.x0(1:(Ne+Ni)*aff,1:(Ne+Ni)*proj); % any (Ne+Ni)*proj rows can be taken in x0
%   STD_AIPtoPMv(find(STD_AIPtoPMv<0)) = 0;
%   
%   % STP PMvto AIP
%   STD_PMvtoAIP = STD_PMvtoAIP +(1-STD_PMvtoAIP)/taud -U*AIP.x0(1:(Ne+Ni)*aff,1:(Ne+Ni)*proj); % any (Ne+Ni)*proj rows can be taken in x0
%   STD_PMvtoAIP(find(STD_PMvtoAIP<0)) = 0;
  
  % plasticity AIP-PMv (any (Ne+Ni)*aff rows can be taken in x0 and x)
  w_PMvtoAIP = w_PMvtoAIP -Ad*PMv.x0(1:(Ne+Ni)*aff,1:(Ne+Ni)*proj).*AIP.v_rekt_md(1+(Ne+Ni)*proj:(Ne+Ni)*proj+(Ne+Ni)*aff) +Ap*PMv.x(1:(Ne+Ni)*aff,1:(Ne+Ni)*proj).*AIP.vy(1+(Ne+Ni)*proj:(Ne+Ni)*proj+(Ne+Ni)*aff).*AIP.v_rekt_mp(1+(Ne+Ni)*proj:(Ne+Ni)*proj+(Ne+Ni)*aff); % weight update
  w_AIPtoPMv = w_AIPtoPMv -Ad*AIP.x0(1:(Ne+Ni)*aff,1:(Ne+Ni)*proj).*PMv.v_rekt_md(1+(Ne+Ni)*proj:(Ne+Ni)*proj+(Ne+Ni)*aff) +Ap*AIP.x(1:(Ne+Ni)*aff,1:(Ne+Ni)*proj).*PMv.vy(1+(Ne+Ni)*proj:(Ne+Ni)*proj+(Ne+Ni)*aff).*PMv.v_rekt_mp(1+(Ne+Ni)*proj:(Ne+Ni)*proj+(Ne+Ni)*aff ); % weight update
  w_PMvtoAIP(w_PMvtoAIP>wmax) = wmax; % weight upper bound
  w_AIPtoPMv(w_AIPtoPMv>wmax) = wmax; % weight upper bound
  
  PMv.x0(:,fired_e_PMv) = 0;
  AIP.x0(:,fired_e_AIP) = 0;
  
  PMv.w(mask_PMv==0) = 0;
  AIP.w(mask_AIP==0) = 0;

%   pot(t) = PMv.v(1);
  wplot(t,:) = [w_PMvtoAIP(1,1) w_PMvtoAIP(1,2) PMv.w(3,2) PMv.w(4,2)];
%   if t==1
%       learn = w_AIPtoPMv;
%   end
%   if t==60e3
%       learn = w_AIPtoPMv-learn;
%   end
end

% figure
% plot(pot);
figure
plot(firings_PMv(:,1),firings_PMv(:,2),'.');
title('PMv')
figure
plot(firings_AIP(:,1),firings_AIP(:,2),'.');
title('AIP')
figure
plot(wplot(:,1));
hold on
plot(wplot(:,2),'r');
plot(wplot(:,3),'g');
plot(wplot(:,4),'k');
hold off

toc;
