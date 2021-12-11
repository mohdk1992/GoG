% Emerging Allee Effect Is Critical For Tumor Persistence/ Benchmark/
% Figure 2 

clear 
clc 

%% Model parameters
L= 10; % Lattice size
r_d = 0.01; % Probabilty of cell dying 
K= 8; % Cell number per node
b = 4; % Number of moving cells channels
r_b = 0.2; % Proliferation probability 
ke = -8:0.1:8; % Switching intensity
Theta = 0:0.05:1; % Switch threshold 
steps = 30; % Time steps

%% Movement
c = [[1,0];[0,1];[-1,0];[0,-1]]; % Moving direction


%% Lattice Initiate

int_radius = 3; % Tumor radious
[xGrid, yGrid] = meshgrid(1:L,1:L);
bw = sqrt((xGrid - (L/2 + 0.5)).^2 + (yGrid - (L/2 + 0.5)).^2) <= int_radius;
bw = double(bw);
Nodes = cat(3,bw,bw); % Initial Lattice

%% Results storage
n_r_int = sum(sum(Nodes(:,:,2)));
n_r_total = zeros(size(Theta));
n_r_k_theta = zeros(length(ke), length(Theta));
y = length(ke);


%% Main
for ke = ke
    e = 1; % Index
    for theta = Theta
        Nodes = cat(3,bw,bw); % Initial Lattice
        for runs = 1:steps  
            % random cell selection
            x = 1:L;
            x_rand = x(randperm(length(x)));
            for i = 1:L
                i = x_rand(i);
                y_rand = x(randperm(length(x)));
                for j = 1:L
                    j = y_rand(j);
                    n = sum(Nodes(i, j, :)); % Total number of cells in the node
                    n_m = Nodes(i, j, 1); % Number of moving cells  in the node
                    n_r = Nodes(i, j, 2); % Number of resting cells in the node
        
                    %% R1 Death stage
                    n_m_death = sum(randsrc(n_m, 1, [1 0;r_d (1-r_d)])); % Number of moving cell death
                    n_r_death = sum(randsrc(n_r, 1, [1 0;r_d (1-r_d)])); % Number of resting cell death
            
                    n_m = n_m - n_m_death;
                    n_r = n_r - n_r_death;
            
                    Nodes(i, j, 1) = n_m;
                    Nodes(i, j, 2) = n_r;
        
                    %% R2 Proliferation of resting cells
                    if n_r > 0
                        for k = 1:n_r
                            proliferation = randsrc(1, 1, [1 0;r_b (1-r_b)]);
                            if n_r < (K - b) && proliferation == 1
                                n_r = n_r + 1;
                            end
                        end
                    end
            
                    Nodes(i, j, 1) = n_m;
                    Nodes(i, j, 2) = n_r;
                    
                    %% R3 Chaning phenotype moving cells
                    if n_m > 0
                        for k = 1:n_m
                            d = (n_m + n_r) / K; % density
                            r_s = 0.5 * (1 + tanh(ke *(d - theta)));
                            sw = randsrc(1, 1, [1 0;r_s (1-r_s)]);
                            if sw == 0 && ke <= 0
                                n_m = n_m - 1;
                                if n_r <(K - b)
                                    n_r = n_r + 1;
                                end
                            end
        
                            if sw == 1 && ke > 0
                                n_m = n_m - 1;
                                if n_r <(K - b)
                                    n_r = n_r + 1;
                                end
                            end
                        end
                    end
                    Nodes(i, j, 1) = n_m;
                    Nodes(i, j, 2) = n_r;

                    %% R3 Chaning phenotype resting cells
                    if n_r > 0
                        for k = 1:n_r
                            d = (n_m + n_r) / K; % density
                            r_s = 0.5 * (1 + tanh(ke *(d - theta)));
                            sw = randsrc(1, 1, [1 0;r_s (1-r_s)]);
                            if sw == 1 && ke <= 0
                                n_r = n_r - 1;
                                if n_m < b
                                    n_m = n_m + 1;
                                end
                            end
        
                            if sw == 0 && ke > 0
                                n_r = n_r - 1;
                                if n_m < b
                                    n_m = n_m + 1;
                                end
                            end
                            
                        end
                    end
                    Nodes(i, j, 1) = n_m;
                    Nodes(i, j, 2) = n_r;
                   %% R4 Moving cell walk
                   if n_m >0
                       for k = 1:n_m
                           ran = [1 2 3 4];
                           direction = c(randsample(ran,1),:);
                           %%%% Boundry condition
                           if (i + direction(1)) > L || (i + direction(1)) <= 0 ||...
                                   (j + direction(2)) > L || (j + direction(2)) <= 0
                               continue
                           end
                           neighbor_n_m = Nodes((i + direction(1)), (j + direction(2)), 1);
                           if neighbor_n_m < b
                               Nodes((i + direction(1)), (j + direction(2)), 1) = neighbor_n_m + 1;
        
                               Nodes(i, j, 1) = Nodes(i, j, 1) - 1;
                           end
                       end
                   end
                   Nodes(i, j, 1) = n_m;
                   Nodes(i, j, 2) = n_r;  
                end
            end
        end
        n_r_new = sum(sum(Nodes(:,:,2)));
        n_r_rate = (n_r_new - n_r_int) / n_r_int;
        n_r_total(e) = n_r_rate;
        e = e + 1;
    end
    n_r_k_theta(y,:) = n_r_total;
    n_r_total = zeros(size(Theta));
    y = y - 1;
end



%% plot heatmap
fig = figure();
ax = axes(fig);
imagesc(ax, (n_r_k_theta/1000 * 6),'XData', 1/2, 'YData', 1/2)
colormap(jet);
colorbar

ke = 8:-0.1:-8;
theta = 0:0.05:1;
Rsize = size(n_r_k_theta);
set(ax, 'xtick', 1:10:Rsize(2), 'xticklabel', [0,0.5,1], ...
'ytick', 1:40:Rsize(1), 'yticklabel', [8, 4, 0, -4, -8])
xlabel("switch threshold (θ)");
ylabel ("switch intensity (k)");
ax.FontSize = 20;

    