%% Clean the workspace and the command window for a fresh start
clearvars;
clc;
%% Reset parpool
delete(gcp);
%% --- Two-traits --- %%

% Parameter Setup
summary.step = 0.1;
summary.sigW = 2;
summary.SelectionBias = [1, 2]; 
summary.InitialAngle = [atan2(1, 6.4), atan2(1, 3.2), atan2(1, 1.6), atan2(1, 0.8), atan2(1, 0.4), atan2(1, 0.2)];
summary.d_ref_List = {-1*[1; 1], -2*[1; 1]};

summary.N = 10^4; 
summary.U = 10^-3;
summary.R = 1;

summary.generations = 15000;
summary.t_burnin = 10000;

% Run simulation
tic
Data = cell(1,length(summary.InitialAngle));
for idx = 1:2
summary.d_ref = summary.d_ref_List{idx};
summary.W_ref = exp(-(vecnorm(summary.d_ref))^2/2/summary.sigW^2);
    
for i_pos = 1:length(summary.InitialAngle)
numrepeat = 100;   
Trait = [1,1]; % [1,1] if two trait model

[ res ] = StdState_R( numrepeat, i_pos, summary, Trait);
j = (idx-1)*length(summary.InitialAngle) + i_pos;
summary.parameter( j,:) = res.parameter;
summary.var1(j,1:numrepeat) = mean(res.var1,2);
summary.var2(j,1:numrepeat) = mean(res.var2,2);
summary.cov(j,1:numrepeat) = mean(res.Covariance,2);

Data{i_pos} = res;
end
end

toc

fname = sprintf('Two_Traits.mat'); 
save(fname, 'Data', 'summary');

delete(gcp)

%% --- Trait 1 only --- %%

% Parameter Setup
summary.step = 0.1;
summary.sigW = 2;
summary.SelectionBias = [1, 2]; 
summary.InitialAngle = [atan2(1, 6.4), atan2(1, 3.2), atan2(1, 1.6), atan2(1, 0.8), atan2(1, 0.4), atan2(1, 0.2)];
summary.d_ref_List = {-1*[1; 1], -2*[1; 1]};

summary.N = 10^4; 
summary.U = 10^-3;
summary.R = 0;

summary.generations = 15000;
summary.t_burnin = 10000;

% Run simulation
tic
Data = cell(1,length(summary.InitialAngle));
for idx = 1:2
summary.d_ref = summary.d_ref_List{idx};
summary.W_ref = exp(-(vecnorm(summary.d_ref))^2/2/summary.sigW^2);
    
for i_pos = 1:length(summary.InitialAngle)
numrepeat = 100;   
Trait = [1,0]; % [1,0] if trait 1 only

[ res ] = StdState_R( numrepeat, i_pos, summary, Trait);
j = (idx-1)*length(summary.InitialAngle) + i_pos;
summary.parameter( j,:) = res.parameter;
summary.var1(j,1:numrepeat) = mean(res.var1,2);
summary.var2(j,1:numrepeat) = mean(res.var2,2);
summary.cov(j,1:numrepeat) = mean(res.Covariance,2);

Data{i_pos} = res;
end
end
toc

fname = sprintf('Trait1.mat'); 
save(fname, 'Data', 'summary');

delete(gcp)

%% --- Trait 2 only --- %%

% Parameter Setup
summary.step = 0.1;
summary.sigW = 2;
summary.SelectionBias = [1, 2]; 
summary.InitialAngle = [atan2(1, 6.4), atan2(1, 3.2), atan2(1, 1.6), atan2(1, 0.8), atan2(1, 0.4), atan2(1, 0.2)];
summary.d_ref_List = {-1*[1; 1], -2*[1; 1]};

summary.N = 10^4; 
summary.U = 10^-3;
summary.R = 0;

summary.generations = 15000;
summary.t_burnin = 10000;

% Run simulation
tic
Data = cell(1,length(summary.InitialAngle));
for idx = 1:2
summary.d_ref = summary.d_ref_List{idx};
summary.W_ref = exp(-(vecnorm(summary.d_ref))^2/2/summary.sigW^2);
    
for i_pos = 1:length(summary.InitialAngle)
numrepeat = 100;   
Trait = [0,1]; % [0,1] if trait 1 only

[ res ] = StdState_R( numrepeat, i_pos, summary, Trait);
j = (idx-1)*length(summary.InitialAngle) + i_pos;
summary.parameter( j,:) = res.parameter;
summary.var1(j,1:numrepeat) = mean(res.var1,2);
summary.var2(j,1:numrepeat) = mean(res.var2,2);
summary.cov(j,1:numrepeat) = mean(res.Covariance,2);

Data{i_pos} = res;
end
end
toc

fname = sprintf('Trait2.mat'); 
save(fname, 'Data', 'summary');

%% Functions

function [ res ] = StdState_R( numrepeat, i_pos, summary, Trait)
%% Set up
k = summary.step;
sigW = summary.sigW;

% Due to unequal selection pressures on modules, the fitness isocline is
% elliptical; here, we find the distance, d, which with an input angle,
% theta, yields the prescribed value of fitness, summary.W_ref.
theta = summary.InitialAngle(i_pos); 
d = Find_InitialPoint( theta, summary.SelectionBias, summary.W_ref, sigW );
d = d(d>0);
WT = d.*[-cos(theta); -sin(theta)];

% rho: Scaling factor for the basal mutation rate
rho = abs(summary.d_ref(1)) + abs(summary.d_ref(2));

a = summary.SelectionBias(1);
b = summary.SelectionBias(2);

N = summary.N;
U = summary.U;
R = summary.R;

%% Enforce the steady state 

F0 = -sqrt(a*WT(1)^2 + b*WT(2)^2)^2/2/sigW^2;
if Trait(1) == 1 && Trait(2) == 1
    % Both traits active
    s1 = (-sqrt(a*(WT(1) + k)^2 + b*WT(2)^2).^2/2/sigW^2) - F0;
    s2 = (-sqrt(a*WT(1)^2 + b*(WT(2) + k)^2).^2/2/sigW^2) - F0;
    U1 = (abs(WT(1))/rho)*U;
    U2 = (abs(WT(2))/rho)*U;
elseif Trait(1) == 1 && Trait(2) == 0
    % Only the first trait is active
    s1 = (-sqrt(a*(WT(1) + k)^2 + b*WT(2)^2).^2/2/sigW^2) - F0;
    s2 = 0;
    U1 = (abs(WT(1))/rho)*U;
    U2 = 0;
elseif Trait(1) == 0 && Trait(2) == 1
    % Only the second trait is active
    s1 = 0;
    s2 = (-sqrt(a*WT(1)^2 + b*(WT(2) + k)^2).^2/2/sigW^2) - F0;
    U1 = 0;
    U2 = (abs(WT(2))/rho)*U;
else
    % Null case
    s1 = 0;
    s2 = 0;
    U1 = 0;
    U2 = 0;
end

%% << Main Loop >>
generations = summary.generations;
t_burnin = summary.t_burnin;

var1 = zeros(numrepeat,generations);
var2 = zeros(numrepeat,generations);
Covariance = zeros(numrepeat,generations);

parfor i_rep = 1:numrepeat
%% --- WF Simulation --- %%
% Population matrix storing meta-information of all individuals.
Popmat = zeros(6,N);

WT = d.*[-cos(theta); -sin(theta)];
Popmat(1,:) = WT(1); % Module 1
Popmat(2,:) = WT(2); % Module 2
Popmat(3,:) = 0; % Number of mutations accumulated in Module 1
Popmat(4,:) = 0; % Number of mutations accumulated in Module 2
Popmat(5,:) = U1; % Mutation rate Module 1
Popmat(6,:) = U2; % Mutation rate Module 2

%% << Main Loop >>

for t = 1:generations
    %% Mutation
    if U1 ~= 0
        nMut_X = poissrnd(sum(Popmat(5,:))); % The number of mutation events
        MutID_X = datasample(1:N, nMut_X, 'Weights', Popmat(5,:),'Replace',false);
        Mutmat_X = Popmat(1:4,MutID_X);
        Mutmat_X(1,:) = Mutmat_X(1,:) + k;
        Mutmat_X(3,:) = Mutmat_X(3,:) + 1;
    
        Popmat(1:4,MutID_X) = Mutmat_X;
    end
    
    if U2 ~= 0
        nMut_Y = poissrnd(sum(Popmat(6,:))); % The number of mutation events
        MutID_Y = datasample(1:N, nMut_Y, 'Weights', Popmat(6,:),'Replace',false);
        Mutmat_Y = Popmat(1:4,MutID_Y);
        Mutmat_Y(2,:) = Mutmat_Y(2,:) + k;
        Mutmat_Y(4,:) = Mutmat_Y(4,:) + 1;
    
        Popmat(1:4,MutID_Y) = Mutmat_Y;
    end

    %% Recombination
    if R ~= 0
        if R == 1  % Full recombination
            nPair = N/2;
        else
            nPair = mybinornd(N, R/2); % The number of recombining pair.
            if 2*nPair > N 
                nPair = N/2;
            end
        end
        
        RecID = randperm(N,2*nPair); % Index of recombining individual.
        Recmat = Popmat(1:4, RecID);
        Recmat(2, 1:nPair) = Recmat(2, nPair+1:2*nPair);
        Recmat(2, nPair+1:2*nPair) = Recmat(2, 1:nPair);
        
        Recmat(4, 1:nPair) = Recmat(4, nPair+1:2*nPair);
        Recmat(4, nPair+1:2*nPair) = Recmat(4, 1:nPair);
        
        Popmat(1:4,RecID) = Recmat;
    end
    
    %% Selection 
    F = Popmat(3,:)*s1 + Popmat(4,:)*s2;
    min_fitness = min(F);
    w = exp( F - min_fitness );
    PROB = w/sum(w);
    
    Smat = zeros(2, N);
    Smat(1,1:N) = 1:N;
    Smat(2,:) = mnrnd(N, PROB);
    
    % Introduce a variable to index number of offsprings.
    Indict = 0; 
    for i = 1:max(Smat(2,:))
        Offsprings = Smat(:, Smat(2,:) == i);
    
        Template = Popmat(:, Offsprings(1,:));
        Template = repmat(Template, 1, i);
        
        Popmat(:, Indict+1 : Indict+length(Template(1,:))) = Template;
        
        Indict = Indict+length(Template(1,:));
    end
    
    %% Meta-information  
    CovMat = cov(Popmat(1,:),Popmat(2,:));
    var1(i_rep,t) = CovMat(1,1);
    var2(i_rep,t) = CovMat(2,2);
    Covariance(i_rep,t) = CovMat(2,1);
    
end
end
res.parameter = [s1;s2;U1;U2;N;R];
res.var1 = var1(:,t_burnin+1:generations);
res.var2 = var2(:,t_burnin+1:generations);
res.Covariance = Covariance(:,t_burnin+1:generations);
end


function d = Find_InitialPoint( theta, ShapeParameter, Default_fitness, sigW )

if nargin < 3
    Default_fitness = 0.000123409804086679;
    sigW = 1;
end

a1 = ShapeParameter(1);
a2 = ShapeParameter(2); 

syms r
eqn = Default_fitness == exp(-sqrt(a1*(r*cos(theta))^2 + a2*(r*sin(theta))^2)^2/2/sigW^2); 
sol = solve(eqn, r);
d = double(sol);
end

% binomial random draw function optimized for speed.
function [ res ] = mybinornd( N, p )
[row_cnt, col_cnt] = size(N);
res = zeros(row_cnt, col_cnt);
for ii=1:row_cnt
   for jj=1:col_cnt
      res(ii, jj) = sum(rand(1,N(ii,jj))<p(ii,jj));
   end
end
end
