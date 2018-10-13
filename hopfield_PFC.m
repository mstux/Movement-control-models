clear

%% Initialization

% number of neurons and actions
N_F5 = 100;
N_AIP = 100;
N_PFC = 100;
N_patterns = 10;

% initialize weights
W_AIP = zeros(N_AIP);
W_rAIPandF5 = zeros(N_AIP/2+N_F5);
W_lAIPandPFC = zeros(N_AIP/2+N_PFC);

% initialize neural patterns corresponding to actions
patterns_F5 = 2*round(rand(N_F5,N_patterns))-1;
patterns_AIP_left = 2*round(rand(N_AIP/2,N_patterns))-1;
patterns_AIP_right = 2*round(rand(N_AIP/2,N_patterns))-1;
patterns_PFC = 2*round(rand(N_PFC,N_patterns))-1;

iter = 10; % number of iterations
similarity = []; % measure of similarity between activity in F5 and motor patterns stored in F5
similarity_test = [];

%% Weights
% Learn synaptic weights using Hopfield rule

patterns_AIP = [patterns_AIP_left; patterns_AIP_right];
for i=1:N_patterns
    W_AIP = W_AIP + patterns_AIP(:,i)*patterns_AIP(:,i)';
end
W_AIP = W_AIP/N_patterns;
W_AIP = W_AIP - diag(diag(W_AIP));
W_AIP(1:N_AIP/2,1:N_AIP/2) = 0;
W_AIP(N_AIP/2+1:end,N_AIP/2+1:end) = 0;

patterns_rAIPandF5 = [patterns_AIP_right; patterns_F5];
for i=1:N_patterns
    W_rAIPandF5 = W_rAIPandF5 + patterns_rAIPandF5(:,i)*patterns_rAIPandF5(:,i)';
end
W_rAIPandF5 = W_rAIPandF5/N_patterns;
W_rAIPandF5 = W_rAIPandF5 - diag(diag(W_rAIPandF5));

patterns_lAIPandPFC = [patterns_AIP_left; patterns_PFC];
for i=1:N_patterns
    W_lAIPandPFC = W_lAIPandPFC + patterns_lAIPandPFC(:,i)*patterns_lAIPandPFC(:,i)';
end
W_lAIPandPFC = W_lAIPandPFC/N_patterns;
W_lAIPandPFC = W_lAIPandPFC - diag(diag(W_lAIPandPFC));
W_lAIPandPFC(N_AIP/2+1:end,1:N_AIP/2) = 0;

%% Dynamics

s = randi(N_patterns); % select pattern
state_PFC = [patterns_PFC(1:N_PFC/2,s); zeros(N_PFC/2,1)]; % half of original pattern;
state_AIP_left = zeros(N_AIP/2,1);
state_AIP_right = zeros(N_AIP/2,1); % 2*round(rand(N_AIP/2,1))-1;
state_F5 = zeros(N_F5,1); % 2*round(rand(N_F5,1))-1;
state = [state_PFC; state_AIP_left; state_AIP_right; state_F5];

for j=1:iter
    % compute activations
    activations_PFC = W_lAIPandPFC(N_AIP/2+1:end,:)*[state_AIP_left; state_PFC];
    activations_AIP_left = W_AIP(1:N_AIP/2,:)*[state_AIP_left; state_AIP_right] + W_lAIPandPFC(1:N_AIP/2,:)*[state_AIP_left; state_PFC];
    activations_AIP_right = W_AIP(N_AIP/2+1:end,:)*[state_AIP_left; state_AIP_right] + W_rAIPandF5(1:N_AIP/2,:)*[state_AIP_right; state_F5];
    activations_F5 = W_rAIPandF5(N_AIP/2+1:end,:)*[state_AIP_right; state_F5];
    activations = [activations_PFC; activations_AIP_left; activations_AIP_right; activations_F5];
    % update
    for i=1:length(activations)
        if activations(i) > 0
            state(i) = 1;
        elseif activations(i) < 0
            state(i) = -1;
        end
    end
    state_PFC =  state(1:N_PFC);
    state_AIP_left = state(N_PFC+1:N_PFC+N_AIP/2);
    state_AIP_right = state(N_PFC+N_AIP/2+1:N_PFC+N_AIP);
    state_F5 = state(N_PFC+N_AIP+1:end);
    
    similarity(end+1) = (state_F5'*patterns_F5(:,s))/N_F5;
    for i=1:N_patterns
        similarity_test(i,j) = (state_F5'*patterns_F5(:,i))/N_F5;
    end
end

figure
hold on
for i=1:N_patterns
    plot(similarity_test(i,:),'r-o')
end
plot(similarity,'b-o')
hold off
