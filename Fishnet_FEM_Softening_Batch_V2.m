
% Unit System: length - mm;     stress - MPa

clear
%% parameter initialization
% dimension
m = 16;                              % number of nacre rows
n = 16;                              % number of nacre colomns 
num_e = m * n;                      % number of elements
num_n = (n+1) * m/2 + n/2;          % number of nodes
connect = conn(m,n);                % connectivity matrix
coord = coords(m,n);                % nodal coordinate vector
b_n = b_node(m,n);                  % boundary node number matrix
in_n = inner_node(b_n,connect,m,n); % internal node number vector
n_1 = numel(b_n);                   % number of Dirichlet B.C nodes
n_2 = num_n - n_1;                  % number of Neuman B.C nodes
start_run = 1;                      % run number of starting simulation
end_run = 1;                       % run number of finishing simulation
num_run = end_run - start_run + 1;  % total number of runs

% geometry
area = 1;                      % cross section area (mm^2) (May be too large, just for simplicity)
length = 10e-3;                % length of element (mm)

% material
elas = 1;                          % Young's modulus (MPa)
n_damage = 1;                      % total number of damages for each link (default = 20)
k0 = elas * area / length;         % initial (original) stiffness of each link
k_t = - 10.0 * k0;                    % tangential stiffness during unloading

% exectution time measurement START
tic;

% nominal strength for each SIMULATION
nom_strength = zeros(num_run,1);

% cluster size for each SIMULATION
s_cluster = zeros(num_run,1);

% gamma-s_cluster
gamma = zeros(num_run,2);

% stress of the Last softened link at the peak step for each SIMULATION
max_sig = zeros(num_run,1);

%% BATCH LOOP

% waitbar
h = waitbar(0,'Percent Completed:');

% Initialization of stiffness matrix K
K_0 = stiff_New(m,n,connect,elas,area,length);


for kk = 1:num_run
    %% variable initialization
    % loading
    u_e_load = 1e-4;               % elongation of a single element
    u_load = u_e_load * n;         % elongation of the whole structure

    %u = zeros(num_n,1);                      % nodal displacement
    disp = zeros(n_damage*num_n+1,1);

    % Initialize u_1
    u_1 = zeros(n_1,1);                      % Dirichlet nodal disp.
    u_1(n_1/2+1:end) ...
            = u_load * ones(n_1/2,1);        
    
    % nodal forces
    f = zeros(num_n,1);                        
    % nominal stress for each INCREMENT
    nominal_sig = zeros(n_damage*num_n+1,1);   
    
    % element stiffness history
    k_vector = elas * area / length * ones(num_e,1);
    % link status
    status_lk = n_damage * ones(num_e,1);
    % link status history (run num. = 1 ONLY)
    if num_run == 1
        status_lk_his = zeros(num_e,num_e);
        sig_his = zeros(num_e,num_e);
        sig_max_his = zeros(2,num_e); 
    end
    
    %% random strength vector
%    sig_rand = r_sig_W_thick(num_e);      % Weibull modulus 5, scale para 10
%    sig_rand = r_sig_G(num_e,10,1);      % Gaussian N(10,2)
%    sig_rand = r_sig_W(num_e);            % Weibull
    sig_rand = r_sig_WG(num_e);             % Weibull-Gaussian Grafted
%    sig_rand = r_sig_PG(num_e);             % Power law - Gaussian

%sig_rand0 = r_sig_WG(num_e);
%sig_rand = sig_rand0(257:512);

srt_r = sort(sig_rand); 

    %% STEP LOOP
    for ii = 1: n_damage*num_e %num_e
        if ii == 1
           K = K_0;
           sig_strength = sig_rand;
        end
        % solve linear system
        [sig,f_1,u] = kernel_New(K,k_vector,connect,area,length,b_n,in_n,n_1,n_2,num_n,num_e,u_1);

        % ratio of random strength over calculated stress
        ratio_sig = sig_strength./sig;
        
        % break condition
        if min(abs(ratio_sig)) > 1e3
            disp(ii+1) = u_1(end); 
            break; 
        end
        % modification of ratio_sig
        j_end = size(ratio_sig);
        j_end = j_end(1);
        for jj = 1:j_end
           if ratio_sig(jj) <= 0
              ratio_sig(jj) = 1e5; 
           end
        end

        % load multiplier
        [k_load,i_fail] = min(ratio_sig);
        % adjust displacement loading
        u_1 = k_load * u_1;
        % solve linear system again
%         [sig,f_1,u] = kernel_New(K,k_vector,connect,area,length,b_n,in_n,n_1,n_2,num_n,num_e,u_1);
%         % check failure
%         ratio_sig = sig_strength./sig;
%         remainder = abs(ratio_sig - 1) < 1e-4;
%         i_fail = find(remainder == 1);
%         i_fail = i_fail(1,1);
        % update true stress and force
        sig = k_load * sig;
        f_1 = k_load * f_1;
        
        % update failure status
        status_lk(i_fail) = status_lk(i_fail) - 1; 
        
        % ************************************************
        % update failure status HISTORY (run num. = 1 ONLY)
%         if num_run == 1
%             if mod(ii,10) == 1
%             sig_his(:,(ii-1)/10+1) = sig;
%             end
%         end
        sig_his(:,ii) = sig;
        %*************************************************
       
        
        % update load and displacement
        tot_Rforce = sum(f_1(1:n_1/2));
        nominal_sig(ii+1) = - tot_Rforce / (area * n_1);
        disp(ii+1) = u_1(end);
        % update strength vector
        sig_strength(i_fail) = (status_lk(i_fail)/n_damage) * sig_rand(i_fail);
        % update stiffness matrix and k_vector
        stiff_old = k_vector(i_fail);
        sig_r = sig_strength(i_fail);
        sig_ori = sig_rand(i_fail);
        stiff_new = stiff_residual(k0,sig_ori,sig_r,k_t);
        K = stiff_Update_s(i_fail,K,connect,stiff_old,stiff_new);
        k_vector(i_fail) = stiff_new;
        
        status_lk_his(:,ii) = status_lk;
        
        
        
        if num_run == 1
           size_damage = sum(status_lk_his(:,ii) ~= n_damage);
           if ii > 1 && size_damage == sum(status_lk_his(:,ii-1) ~= n_damage)
               if (- tot_Rforce / (area * n_1)) > temp
                   sig_max_his(1,size_damage) = - tot_Rforce / (area * n_1);
                   sig_max_his(2,size_damage) = srt_r(size_damage);%*(num_e-size_damage)/num_e;
                   temp = - tot_Rforce / (area * n_1);
               end
               
           else
               sig_max_his(1,size_damage) = - tot_Rforce / (area * n_1);
               sig_max_his(2,size_damage) = srt_r(size_damage);%*(num_e-size_damage)/num_e;
               temp = - tot_Rforce / (area * n_1);
           end
        end
        
%         if num_run == 1 && ii <= num_e
%              sig_max_his(1,ii) = - tot_Rforce / (area * n_1);
%              sig_max_his(2,ii) = srt_r(ii);%*(num_e-size_damage)/num_e; 
%         end
        % data needed by the figure in the paper
%         if ii == 173
%             size_damage
%             i_fail
%         end
        
    end % end of step loop

    
    ind_max = find(nominal_sig == max(nominal_sig));
    % peak nominal stress
    nom_strength(kk) = max(nominal_sig);
    max_sig(kk,1) = max(sig_his(:,ind_max));
    
    %plot(disp,nominal_sig);
    %figure;
    %plot(nominal_sig./disp);
    %figure;
    %plot(disp);
    %figure;
    %plot(nominal_sig);

     % damage zone size
     s_cluster(kk) = sum(status_lk_his(:,ind_max) ~= n_damage);


%     plot load-displacement curve
    if num_run == 1
        figure
        plot(disp,nominal_sig)
        figure
 %       plot up to maximum load
       plot(sig_max_his(2,1:s_cluster)./sig_max_his(1,1:s_cluster),'--')
 %       plot up to the end of localization
       plot(sig_max_his(2,:)./sig_max_his(1,:))
       hold on;
    end
    
%     if num_run == 100
%           gamma(kk,1) = s_cluster(kk);
%           gamma(kk,2) = sig_max_his(2,s_cluster(kk))./sig_max_his(1,s_cluster(kk)); 
%     end
%     plot(sig_max_his(1,1:s_cluster(kk)))
%     hold on;
%     plot(sig_max_his(2,1:s_cluster(kk)))

    
    % Progress
    progress = kk / num_run;
    waitbar(progress)
end % end of batch loop
close(h)
% exectution time measurement END
toc;
%%
%BATCH output

% if num_run > 1
%     filename = '------FILE PATH-------';
%     sheet = 'Sheet1';
%     xlRange = 'A1';
%     %xlswrite(filename,[disp,nominal_sig],sheet,xlRange)
%     xlswrite(filename,nom_strength,sheet,xlRange);
%     %save('Nacre_SE.mat','nominal_sig')
%     
%   
%     filename = '------FILE PATH-------';
%     sheet = 'Sheet1';
%     xlRange = 'A1';
%     xlswrite(filename,s_cluster,1,xlRange)
% 
% end


%% Single Realization Output

% *************************************
%   Make sure run number is set to 1
% *************************************

% if num_run == 1
%     filename = '------FILE PATH-------';
%     sheet = 'Sheet1';
%     xlRange = 'A1';
%     xlswrite(filename,coord,1,xlRange)
% 
%     filename = '------FILE PATH-------';
%     sheet = 'Sheet1';
%     xlRange = 'A1';
%     xlswrite(filename,connect,1,xlRange)
% 
%     filename = '------FILE PATH-------';
%     sheet = 'Sheet1';
%     xlRange = 'A1';
%     xlswrite(filename,status_lk_his,1,xlRange)
% 
%     filename = '------FILE PATH-------';
%     sheet = 'Sheet1';
%     xlRange = 'A1';
%     xlswrite(filename,sig_his,1,xlRange)
%     % 
%     filename = '------FILE PATH-------';
%     sheet = 'Sheet1';
%     xlRange = 'A1';
%     xlswrite(filename,sig_rand,1,xlRange)
%
% end
%%
%     filename = '------FILE PATH-------';
%     sheet = 'Sheet1';
%     xlRange = 'A1';
%     xlswrite(filename,coord,1,xlRange)
% 
%     filename = '------FILE PATH-------';
%     sheet = 'Sheet1';
%     xlRange = 'A1';
%     xlswrite(filename,connect,1,xlRange)
% 
%     filename = '------FILE PATH-------';
%     sheet = 'Sheet1';
%     xlRange = 'A1';
%     xlswrite(filename,status_lk_his,1,xlRange)
% 
%     filename = '------FILE PATH-------';
%     sheet = 'Sheet1';
%     xlRange = 'A1';
%     xlswrite(filename,sig_his,1,xlRange)
%     % 
%     filename = '------FILE PATH-------';
%     sheet = 'Sheet1';
%     xlRange = 'A1';
%     xlswrite(filename,sig_rand,1,xlRange)
% 
%     filename = '------FILE PATH-------';
%     sheet = 'Sheet1';
%     xlRange = 'A1';
%     xlswrite(filename,[disp,nominal_sig],sheet,xlRange)

%% post-processing (visualization)
% stress level
% for jj = 1:ii%find(nominal_sig==max(nominal_sig))
%     if mod(jj,10) == 0
%         plt_ele_scalar(sig_his,num_e,connect,coord,m,n,jj);
%         pause(0.0001);
%         jj/ii
%     end
% end

% damage level
% for jj = 1:find(nominal_sig==max(nominal_sig))
% for jj = 1: ii-1
%     if mod(jj,10) == 0
%         plt_ele_scalar(n_damage-status_lk_his,num_e,connect,coord,m,n,jj);
%         pause(0.0001);
%         jj/ii
%     end
% end
% 
% mean(s_cluster)