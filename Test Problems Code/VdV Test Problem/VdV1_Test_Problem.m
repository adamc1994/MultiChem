% VdV1 test problem.

%% Prompt user to select algorithm file

            [file,path] = uigetfile('*.m');
            selectedfile = fullfile(path,file);

%% Generate initial dataset

% Latin Hypercube (LHC)

        %LHC prep
        
            initial_sample_size = 20;
            inputdim = 2;
            initalX = lhsdesign(initial_sample_size,inputdim);
            lowerbounds = [0.5 25]; % time, temp
            upperbounds = [10 100]; % time, temp          
            
            for i = 1 : inputdim               
                
                initalX(:,i) = initalX(:,i) .*  (upperbounds(i) - lowerbounds(i) ) + lowerbounds(i);
                
            end
            
            conditions = initalX;
            tR = conditions(:,1); % in min 
            T = conditions(:,2); % in oC, rate function converts
            ConcA = 1;
            
        % end LHC prep  
            
%% Calculate converions for LHC experiments

            SM_yield = zeros((size(conditions(:,1),1)) , 1);
            product_yield = zeros((size(conditions(:,1),1)) , 1);
            side1_yield = zeros((size(conditions(:,1),1)) , 1);
            side2_yield = zeros((size(conditions(:,1),1)) , 1);
            E_factor = zeros((size(conditions(:,1),1)) , 1);
            STY = zeros((size(conditions(:,1),1)) , 1);
            
for i=1:1:(size(conditions(:,1),1))
    
            [Cfinal, Yield, TimeData, ConcData] = VdV_Constants(ConcA, T(i) * ones(1,4), tR(i));
            %lastEntry = size(TimeData);
            %lastEntry = lastEntry(1,1);

            % Calculate percentage yields
            
            SM = Yield(1)*100;
            product = Yield(2)*100;
            side1 = Yield(3)*100;
            side2 = Yield(4)*100;
            
            % Adjust yields with error
            
            SM = SM+((rand-0.5)*0.5)+(SM*((rand-0.5)/100));
            if SM < 0
                SM = 1*10^-6;
            end
            if SM > 100
                SM = 100;
            end
            
            product = product+((rand-0.5)*0.5)+(product*((rand-0.5)/100));
            if product < 0
                product = 1*10^-6;
            end
            if product > 100
                product = 100;
            end
            
            side1 = side1+((rand-0.5)*0.5)+(side1*((rand-0.5)/100));
            if side1 < 0
                side1 = 1*10^-6;
            end
            if side1 > 100
                side1 = 100;
            end
            
            side2 = side2+((rand-0.5)*0.5)+(side2*((rand-0.5)/100));
            if side2 < 0
                side2 = 1*10^-6;
            end
            if side2 > 100
                side2 = 100;
            end   
            
            % Calculate other objectives
            
            STY_ = STY_calc(product/100, tR(i), ConcData(1,1), 100);
            
            % Update data
            
            SM_yield(i) = SM;
            product_yield(i) = product;
            side1_yield(i) = side1;
            side2_yield(i) = side2;
            STY(i) = STY_;
                     
            f = [-log(product_yield) -log(STY)]; % maximise yield, maximise STY
                
% Update data set

            data.x(i,:) = conditions(i,:);
            data.y(i,:) = f(i,:);   

% Determine number of Pareto optimal solutions after each iteration
            
            Product_yield = exp(-data.y(:,1));
            STY_actual = exp(-data.y(:,2));

            data.z = [Product_yield STY_actual];
            [opt, data.idxs] = paretoFront(data.z);
            data.optconds = data.x(data.idxs,:); % extract conditions corresponding to Pareto optimal solutions
            opt_size = size(opt);
            opt_n = opt_size(1,1);
            data.opt(i,:) = opt_n;
            data.front = opt;
            
% Calculate hypervolume after each iteration (requires known reference points)

            F = opt; 
            AU = [3.8199 159.9906];  
            U = [31.3995 560.6544];
            N = 100000;
            hv = hypervolume(F, AU, U, N).*100;
            data.hv(i,:) = hv;
            
% Plot hypervolume

            plot(1:i, data.hv, 'LineWidth', 2)
            title('Hypervolume')
            xlabel('Experiment Number')
            ylabel('Hypervolume')
            
            pause(1)
            
end
            
%% Run algorithms
           
% Termination criteria (after xx experiments)

            total_expts = (size(data.x,1));
            
            while total_expts < 100   
                             
% Run desired algorithm

            run(selectedfile);
            
            tR = conditions(:,1); % in min as pre-exp factors (A) are in min-1
            T = conditions(:,2); % in oC, rate function converts
            ConcA = 1;
            
            SM_yield = zeros((size(conditions(:,1),1)) , 1);
            product_yield = zeros((size(conditions(:,1),1)) , 1);
            side1_yield = zeros((size(conditions(:,1),1)) , 1);
            side2_yield = zeros((size(conditions(:,1),1)) , 1);
            E_factor = zeros((size(conditions(:,1),1)) , 1);
            STY = zeros((size(conditions(:,1),1)) , 1);
            
for i=1:1:(size(conditions(:,1),1))
    
           [Cfinal, Yield, TimeData, ConcData] = VdV_Constants(ConcA, T(i) * ones(1,4), tR(i));
            
            % Calculate percentage yields
            
            SM = Yield(1)*100;
            product = Yield(2)*100;
            side1 = Yield(3)*100;
            side2 = Yield(4)*100;
            
            % Adjust yields with error
            
            SM = SM+((rand-0.5)*0.5)+(SM*((rand-0.5)/100));
            if SM < 0
                SM = 1*10^-6;
            end
            if SM > 100
                SM = 100;
            end
            
            product = product+((rand-0.5)*0.5)+(product*((rand-0.5)/100));
            if product < 0
                product = 1*10^-6;
            end
            if product > 100
                product = 100;
            end
            
            side1 = side1+((rand-0.5)*0.5)+(side1*((rand-0.5)/100));
            if side1 < 0
                side1 = 1*10^-6;
            end
            if side1 > 100
                side1 = 100;
            end
            
            side2 = side2+((rand-0.5)*0.5)+(side2*((rand-0.5)/100));
            if side2 < 0
                side2 = 1*10^-6;
            end
            if side2 > 100
                side2 = 100;
            end   
            
            % Calculate other objectives
            
            STY_ = STY_calc(product/100, tR(i), ConcData(1,1), 100);
            
            % Update data
            
            SM_yield(i) = SM;
            product_yield(i) = product;
            side1_yield(i) = side1;
            side2_yield(i) = side2;
            STY(i) = STY_;
            
            f = [-log(product_yield) -log(STY)]; % maximise yield, maximise STY
            
            % Update data set
                
            data.x = [data.x ; conditions(i,:)];
            data.y = [data.y ; f(i,:)];
            total_expts = (size(data.x,1)); % calculate total experiments so far
            
            % Determine number of Pareto optimal solutions after each iteration
            
            Product_yield = exp(-data.y(:,1));
            STY_actual = exp(-data.y(:,2));

            data.z = [Product_yield STY_actual];
            [opt, data.idxs] = paretoFront(data.z);
            data.optconds = data.x(data.idxs,:); % extract conditions corresponding to Pareto optimal solutions
            opt_size = size(opt);
            opt_n = opt_size(1,1);
            data.opt = [data.opt ; opt_n];
            data.front = opt;
            
            % Calculate hypervolume after each iteration (requires known reference points)

            F = opt; 
            AU = [3.8199 159.9906];  
            U = [31.3995 560.6544];
            N = 100000;
            hv = hypervolume(F, AU, U, N).*100;
            data.hv = [data.hv ; hv];
            
% Plot hypervolume

            plot(1:total_expts, data.hv, 'LineWidth', 2)
            title('Hypervolume')
            xlabel('Experiment Number')
            ylabel('Hypervolume')
            
            pause(1)
         
end
            end