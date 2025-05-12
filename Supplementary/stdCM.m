%% Clean the workspace and the command window for a fresh start
clearvars;
clc;
%% Reset parpool
delete(gcp);
%% --- Parameter Setup --- %
% Initialize the population size N
Parameter.N = 10^4;
n_v = 6;

Parameter.v1 = exp(linspace(log(5*10^-6), log(1*10^-2), n_v));
Parameter.v2 = flip(exp(linspace(log(5*10^-6), log(1*10^-2), n_v)))';

% Create the full meshgrid
[X, Y] = meshgrid(Parameter.v1, Parameter.v2);

% Focus on the regime where v1 < v2
M = Y >= X;

X = M.*X;
v1 = X(X > 0);

Y = M.*Y;
v2 = Y(Y > 0);

% Set a, where a is s = a*U. a must be >> 1
n_a = 4;
a = round(exp(linspace(log(1*10^1), log(1*10^2), n_a)));
Parameter.a = a;

% Create Parameter Space
N = Parameter.N;
N0 = N;

tic
for i = 1:numel(v1)
    target_v1 = v1(i);
    target_v2 = v2(i);
    for j = 1:numel(a)
        a0 = a(j);
        syms U0
        f_DF_v1 = ((a0*U0)^2 * (2*log(N0*a0*U0) - log(a0))) / (log(a0))^2 - target_v1;
        f_DF_v2 = ((a0*U0)^2 * (2*log(N0*a0*U0) - log(a0))) / (log(a0))^2 - target_v2;
        
        % Solve the equation numerically
        U1 = double( vpasolve(f_DF_v1, U0, [0, 1]) );
        s1 = a0*U1;
        
        U2 = double( vpasolve(f_DF_v2, U0, [0, 1]) );
        s2 = a0*U2;
        
        % Save parameters in cell array
        Parameter.s1U1{i,j} = [s1;U1];
        Parameter.s2U2{i,j} = [s2;U2];
    end
end

% Check the accuracy
for i = 1:numel(v1)
    syms s0 U0 N0
    AdaptRate = ( (s0^2)*( (2*log(N0*s0))-log(s0/U0) ) )/( log(s0/U0) )^2;
    for j = 1:numel(Parameter.a)
        s1U1 = Parameter.s1U1{i,j};
        s2U2 = Parameter.s2U2{i,j};
        
        Parameter.v1_check(i,j) = double(subs(AdaptRate, {s0, U0, N0}, {s1U1(1), s1U1(2), N}));
        Parameter.v2_check(i,j) = double(subs(AdaptRate, {s0, U0, N0}, {s2U2(1), s2U2(2), N}));
    end
end
toc

%% --- WF Simulation --- %%
tic
parpool;

% Set parameters for the simulation
generations = 15000;
t_burnin = 10000;

[ID_1, ID_2] = meshgrid(1:numel(a), 1:numel(a));

WFsim.generations = generations;

WFsim.t_burnin = t_burnin;

WFsim.F1_2 = cell(numel(v1), 1);
WFsim.F2_1 = cell(numel(v1), 1);

WFsim.v1_2 = cell(numel(v1), 1);
WFsim.v2_1 = cell(numel(v1), 1);

% Set broadcast variables
% Record the whole fitness trajectory for each trait
WFsim_F1_2_all = WFsim.F1_2;
WFsim_F2_1_all = WFsim.F2_1;
% Record the mean rate of adaptation v_{1|2} and v_{2|1}
WFsim_v1_all = WFsim.v1_2;
WFsim_v2_all = WFsim.v2_1;

% Transform 2D cell arrays into 1D cell arrays
v1_param_all = Parameter.s1U1(:,1:end);
v2_param_all = Parameter.s2U2(:,1:end);
v1_param_all_1D = arrayfun(@(i) v1_param_all(i, :), 1:size(v1_param_all, 1), 'UniformOutput', false);
v2_param_all_1D = arrayfun(@(i) v2_param_all(i, :), 1:size(v2_param_all, 1), 'UniformOutput', false);

parfor i = 1:numel(v1)
    % Each parameter set in the same row generates the same value of v1 and
    % v2, respectively.
    v1_param = v1_param_all_1D{i};
    v2_param = v2_param_all_1D{i};
    
    Rec_F1_2 = zeros(numel(ID_1), generations);
    Rec_F2_1 = zeros(numel(ID_2), generations);
    
    Rec_v1  = zeros(1,numel(ID_1));
    Rec_v2  = zeros(1,numel(ID_2));
    
    for j = 1:numel(ID_1)
        % Set s1, U1, s2, and U2
        s1 = v1_param{ID_1(j)}(1);
        U1 = v1_param{ID_1(j)}(2);
        s2 = v2_param{ID_2(j)}(1);
        U2 = v2_param{ID_2(j)}(2);
        
        % Initial population starts with one strain: N individuals with no
        % accumulated mutations in trait 1 and 2.        
        pop = [N 0 0];
        for t = 1:generations
            % Compute the total fitness of each strain  
            min_fitness = min( (pop(:,2)*s1) + (pop(:,3)*s2) );
            w = exp( ((pop(:,2)*s1) + (pop(:,3)*s2)) - min_fitness );
            % Normalize fitness   
            fitness = pop(:,1).* w; 
            PROB = fitness/sum(fitness);
            % Use Poisson sampling to choose parents for next generation    
            parents_counts = poissrnd(N * PROB);
            
            % Mutate population   
            U = U1+U2;
            tot_mutations = binornd(parents_counts, U); 
            num_mutations1 = binornd(tot_mutations, U1/U); 
            parents_counts = parents_counts - num_mutations1;
            num_mutations2 = tot_mutations - num_mutations1; 
            parents_counts = parents_counts - num_mutations2;
            
            mutation1_indices = (num_mutations1 > 0); 
            mutation2_indices = (num_mutations2 > 0);      
            
            pop_new = [parents_counts, pop(:, 2:3)];      
            pop_new1 = zeros(0, 3);
            pop_new2 = zeros(0, 3); 
            
            if any(mutation1_indices) 
               pop_new1 = [num_mutations1(mutation1_indices), pop(mutation1_indices, 2) + 1, pop(mutation1_indices, 3)];       
            end
            if any(mutation2_indices)
               pop_new2 = [num_mutations2(mutation2_indices), pop(mutation2_indices, 2), pop(mutation2_indices, 3) + 1];    
            end
            pop = [pop_new; pop_new1; pop_new2];
            
            % Remove extinct lineages (rows with first column == 0)  
            pop(pop(:, 1) == 0, :) = [];
            
            % Record average fitness per trait       
            Rec_F1_2(j, t) = sum(pop(:, 1) .* pop(:, 2) * s1 ) / sum(pop(:,1));
            Rec_F2_1(j, t) = sum(pop(:, 1) .* pop(:, 3) * s2 ) / sum(pop(:,1));
            
            if t == generations   
                Rec_v1(j) = (Rec_F1_2(j, generations) - Rec_F1_2(j, t_burnin)) / (generations-t_burnin+1); 
                Rec_v2(j) = (Rec_F2_1(j, generations) - Rec_F2_1(j, t_burnin)) / (generations-t_burnin+1);
            end            
        end
    end
    
    WFsim_F1_2_all{i} = Rec_F1_2;
    WFsim_F2_1_all{i} = Rec_F2_1;
    
    WFsim_v1_all{i} = Rec_v1';
    WFsim_v2_all{i} = Rec_v2';
end

WFsim.F1_2 = WFsim_F1_2_all;
WFsim.F2_1 = WFsim_F2_1_all;

WFsim.v1_2 = WFsim_v1_all;
WFsim.v2_1 = WFsim_v2_all;

toc

%% --- WF Simulation, v' estimate --- %%
tic

% Set parameters for the simulation
generations = 15000;
t_burnin = 10000;

WFsim_vi.generations = generations;

WFsim_vi.t_burnin = t_burnin;

vectorStrings = Parameter.s1U1;
v_prime_param_all(1,1:4) = vectorStrings(1,1:4);
v_prime_param_all(2,1:4) = vectorStrings(7,1:4);
v_prime_param_all(3,1:4) = vectorStrings(12,1:4);
v_prime_param_all(4,1:4) = vectorStrings(16,1:4);
v_prime_param_all(5,1:4) = vectorStrings(19,1:4);
v_prime_param_all(6,1:4) = vectorStrings(21,1:4);
v_prime_param_all = v_prime_param_all';

WFsim_vi.F1 = cell(numel(v_prime_param_all), 1);
WFsim_vi.v1 = cell(numel(v_prime_param_all), 1);

% Set broadcast variables
WFsim_F1_all = WFsim_vi.F1;
WFsim_v1_all = WFsim_vi.v1;

parfor i = 1:numel(v_prime_param_all)
    Rec_F1 = zeros(1, generations);
    Rec_v1  = zeros(1);

    % Set s1, U1, s2, and U2
    s1 = v_prime_param_all{i}(1);
    U1 = v_prime_param_all{i}(2);

    % Initial population starts with one strain: N individuals with no
    % accumulated mutations in trait 1 and 2.
    pop = [N 0 0];
    for t = 1:generations        
        % Compute the total fitness of each strain    
        min_fitness = min( (pop(:,2)*s1) );  
        w = exp( (pop(:,2)*s1) - min_fitness );    
            
        % Normalize fitness       
        fitness = pop(:,1).* w;    
        PROB = fitness/sum(fitness);    
        % Use Poisson sampling to choose parents for next generation       
        parents_counts = poissrnd(N * PROB);  
               
        % Mutate population               
        num_mutations1 = binornd(parents_counts, U1);      
        parents_counts = parents_counts - num_mutations1;
        mutation1_indices = (num_mutations1 > 0);                      
        pop_new = [parents_counts, pop(:, 2:3)]; 
        
        % Initialize with correct column size             
        pop_new1 = zeros(0, 3);
        
        if any(mutation1_indices)        
            pop_new1 = [num_mutations1(mutation1_indices), pop(mutation1_indices, 2) + 1, pop(mutation1_indices, 3)];                  
        end        
        pop = [pop_new; pop_new1];
        
        % Remove extinct lineages (rows with first column == 0)  
        pop(pop(:, 1) == 0, :) = [];     
        
        % Record average fitness per trait                 
        Rec_F1(1, t) = sum(pop(:, 1) .* pop(:, 2) * s1 ) /sum(pop(:,1));
           
        if t == generations                  
            Rec_v1 = (Rec_F1(1, generations) - Rec_F1(1, t_burnin)) / (generations-t_burnin+1);            
        end        
    end

    WFsim_F1_all{i} = Rec_F1;
    WFsim_v1_all{i} = Rec_v1';
end

WFsim_vi.F1 = WFsim_F1_all;
WFsim_vi.v1 = WFsim_v1_all;

toc

%% --- Analytical Prediction, equal weights --- %%
% Begin timing the execution of this block
tic
% Define symbolic variables for analytical calculation
syms s0 U0 N0
% Define the mathematical expression for Adaptation rate
AdaptRate = ( (s0^2)*( (2*log(N0*s0))-log(s0/U0) ) )/( log(s0/U0) )^2;

% DF model assumes s >> U, thus, we set tolerance limit "tol" to ensure s
% is at least s = tol*U.
tol = 10;

% Initialize vectors for storing v1_2 and v2_1 calculation results
v1_2 = zeros(numel(WFsim.F1_2),numel(ID_1),1);
v2_1 = zeros(numel(WFsim.F2_1),numel(ID_2),1);

N = Parameter.N;
for i = 1:numel(WFsim.F1_2)
    for j = 1:numel(ID_2)
        v2 = Parameter.v2_check(i,ID_2(j));
        s1 = Parameter.s1U1{i,ID_1(j)}(1);
        U1 = Parameter.s1U1{i,ID_1(j)}(2);
        
        v1 = Parameter.v1_check(i,ID_1(j));
        s2 = Parameter.s2U2{i,ID_2(j)}(1);
        U2 = Parameter.s2U2{i,ID_2(j)}(2);
        
            s_tilde = mean([s1;s2]);
            U2_tilde = double(solve(subs(AdaptRate, {s0, N0}, {s_tilde, N}) == v2, U0));
            U2_tilde = U2_tilde(s_tilde./U2_tilde > tol);
            
            U1_tilde = double(solve(subs(AdaptRate, {s0, N0}, {s_tilde, N}) == v1, U0));
            U1_tilde = U1_tilde(s_tilde./U1_tilde > tol);
            
            if isempty(U1_tilde) == 0 && isempty(U2_tilde) == 0         
                temp_v1_2 = double((U1_tilde/(U1_tilde+U2_tilde))*subs(AdaptRate, {s0, U0, N0}, {s_tilde, U1_tilde+U2_tilde, N}));
                temp_v2_1 = double((U2_tilde/(U1_tilde+U2_tilde))*subs(AdaptRate, {s0, U0, N0}, {s_tilde, U1_tilde+U2_tilde, N}));
            else
                temp_v1_2 = NaN;
                temp_v2_1 = NaN;
            end       
            
        v1_2(i,j) = (temp_v1_2);
        v2_1(i,j) = (temp_v2_1);
    end
end

% Save the results of v1_2 and v2_1 into the Prediction structure
Prediction.v1_2 = v1_2;
Prediction.v2_1 = v2_1;

%% --- Analytical Prediction, weights = U_i/U --- %%
% Begin timing the execution of this block
tic

% Initialize vectors for storing v1_2 and v2_1 calculation results
v1_2 = zeros(numel(WFsim.F1_2),numel(ID_1),1);
v2_1 = zeros(numel(WFsim.F2_1),numel(ID_2),1);

N = Parameter.N;
for i = 1:numel(WFsim.F1_2)
    for j = 1:numel(ID_2)
        v1 = Parameter.v1_check(i,ID_1(j));
        s2 = Parameter.s2U2{i,ID_2(j)}(1);
        U2 = Parameter.s2U2{i,ID_2(j)}(2);
        
        v2 = Parameter.v2_check(i,ID_2(j));
        s1 = Parameter.s1U1{i,ID_1(j)}(1);
        U1 = Parameter.s1U1{i,ID_1(j)}(2);
        
        U_tot = U1 + U2;
        
            s_tilde = (U1/U_tot)*s1+(U2/U_tot)*s2; %mean([s1;s2]);
            U2_tilde = double(solve(subs(AdaptRate, {s0, N0}, {s_tilde, N}) == v2, U0));
            U2_tilde = U2_tilde(s_tilde./U2_tilde > tol);
            
            U1_tilde = double(solve(subs(AdaptRate, {s0, N0}, {s_tilde, N}) == v1, U0));
            U1_tilde = U1_tilde(s_tilde./U1_tilde > tol);
            
            if isempty(U1_tilde) == 0 && isempty(U2_tilde) == 0         
                temp_v1_2 = double((U1_tilde/(U1_tilde+U2_tilde))*subs(AdaptRate, {s0, U0, N0}, {s_tilde, U1_tilde+U2_tilde, N}));
                temp_v2_1 = double((U2_tilde/(U1_tilde+U2_tilde))*subs(AdaptRate, {s0, U0, N0}, {s_tilde, U1_tilde+U2_tilde, N}));
            else
                temp_v1_2 = NaN;
                temp_v2_1 = NaN;
            end       
            
        v1_2(i,j) = (temp_v1_2);
        v2_1(i,j) = (temp_v2_1);
    end
end

% Save the results of v1_2 and v2_1 into the Prediction structure
Prediction_U_weights.v1_2 = v1_2;
Prediction_U_weights.v2_1 = v2_1;

%% --- Analytical Prediction, weights = s_i/(s_1 + s_2) --- %%
tic

% Initialize vectors for storing v1_2 and v2_1 calculation results
v1_2 = zeros(numel(WFsim.F1_2),numel(ID_1),1);
v2_1 = zeros(numel(WFsim.F2_1),numel(ID_2),1);

N = Parameter.N;
for i = 1:numel(WFsim.F1_2)
    for j = 1:numel(ID_2)  
        v1 = Parameter.v1_check(i,ID_1(j));
        s2 = Parameter.s2U2{i,ID_2(j)}(1);
        U2 = Parameter.s2U2{i,ID_2(j)}(2);
        
        v2 = Parameter.v2_check(i,ID_2(j));
        s1 = Parameter.s1U1{i,ID_1(j)}(1);
        U1 = Parameter.s1U1{i,ID_1(j)}(2);
        
        s_tot = s1 + s2;
        
            s_tilde = (s1/s_tot)*s1+(s2/s_tot)*s2;  
            U2_tilde = double(solve(subs(AdaptRate, {s0, N0}, {s_tilde, N}) == v2, U0));
            U2_tilde = U2_tilde(s_tilde./U2_tilde > tol);
            
            U1_tilde = double(solve(subs(AdaptRate, {s0, N0}, {s_tilde, N}) == v1, U0));
            U1_tilde = U1_tilde(s_tilde./U1_tilde > tol);
            
            if isempty(U1_tilde) == 0 && isempty(U2_tilde) == 0         
                temp_v1_2 = double((U1_tilde/(U1_tilde+U2_tilde))*subs(AdaptRate, {s0, U0, N0}, {s_tilde, U1_tilde+U2_tilde, N}));
                temp_v2_1 = double((U2_tilde/(U1_tilde+U2_tilde))*subs(AdaptRate, {s0, U0, N0}, {s_tilde, U1_tilde+U2_tilde, N}));
            else
                temp_v1_2 = NaN;
                temp_v2_1 = NaN;
            end       
            
        v1_2(i,j) = (temp_v1_2);
        v2_1(i,j) = (temp_v2_1);
    end
end

% Save the results of v1_2 and v2_1 into the Prediction structure
Prediction_s_weights.v1_2 = v1_2;
Prediction_s_weights.v2_1 = v2_1;

% End timing and print the elapsed time to Command Window
toc

%% --- Save Data --- %%
fname = sprintf('Stdstate_N_%d.mat', N);
save(fname, 'WFsim','WFsim_vi','Parameter','Prediction', 'Prediction_U_weights','Prediction_s_weights');
