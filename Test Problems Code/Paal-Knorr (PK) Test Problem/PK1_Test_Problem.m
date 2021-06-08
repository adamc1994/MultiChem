% Paal-Knorr1 test problem.

%% Prompt user to select algorithm file

            [file,path] = uigetfile('*.m');
            selectedfile = fullfile(path,file);

%% Generate initial dataset

% Latin Hypercube (LHC)

        %LHC prep
        
            initial_sample_size = 20;
            inputdim = 2;
            initalX = lhsdesign(initial_sample_size,inputdim);
            lowerbounds = [0.5 1]; % time, equiv B
            upperbounds = [2 10]; % time, equiv B
                        
            for i = 1 : inputdim               
                
                initalX(:,i) = initalX(:,i) .*  (upperbounds(i) - lowerbounds(i) ) + lowerbounds(i);
                
            end
            
            conditions = initalX;
            tR = conditions(:,1); % in min 
            T = 50; % in oC, rate function converts
            ConcA = 1;
            equivB = conditions(:,2);
            ConcB = ConcA.*equivB;
            
        % end LHC prep  
            
%% Calculate converions for LHC experiments

            SM_yield = zeros((size(conditions(:,1),1)) , 1);
            SM2_yield = zeros((size(conditions(:,1),1)) , 1);
            intermediate_yield = zeros((size(conditions(:,1),1)) , 1);
            product_yield = zeros((size(conditions(:,1),1)) , 1);
            E_factor = zeros((size(conditions(:,1),1)) , 1);
            STY = zeros((size(conditions(:,1),1)) , 1);
            RME = zeros((size(conditions(:,1),1)) , 1);
            
for i=1:1:(size(conditions(:,1),1))
    
            [Cfinal, Yield, TimeData, ConcData] = PK_Constants(ConcA, ConcB(i), T * ones(1,4), tR(i));
            %lastEntry = size(TimeData);
            %lastEntry = lastEntry(1,1);

            % Calculate percentage yields
            
            SM = Yield(1)*100;
            SM2 = Yield(2)*100;
            intermediate = Yield(3)*100;
            product = Yield(4)*100;
            
            % Adjust yields with error
            
            SM = SM+((rand-0.5)*0.5)+(SM*((rand-0.5)/100));
            if SM < 0
                SM = 1*10^-6;
            end
            if SM > 100
                SM = 100;
            end
            
            SM2 = SM2+((rand-0.5)*0.5)+(SM2*((rand-0.5)/100));
            if SM2 < 0
                SM2 = 1*10^-6;
            end
            if SM2 > 100
                SM2 = 100;
            end
            
            intermediate = intermediate+((rand-0.5)*0.5)+(intermediate*((rand-0.5)/100));
            if intermediate < 0
                intermediate = 1*10^-6;
            end
            if intermediate > 100
                intermediate = 100;
            end 
            
            product = product+((rand-0.5)*0.5)+(product*((rand-0.5)/100));
            if product < 0
                product = 1*10^-6;
            end
            if product > 100
                product = 100;
            end
            
            % Calculate other objectives
            
            STY_ = STY_calc(product/100, tR(i), ConcData(1,1), 139.20);

            RME_ = (139.20*product)/(114.14+(61.08*(ConcData(1,5+1)./ConcData(1,1)))); % 5 = global 'n' (see Lactose constants)
            
            % Update data
            
            SM_yield(i) = SM;
            SM2_yield(i) = SM2;
            intermediate_yield(i) = intermediate;
            product_yield(i) = product;
            STY(i) = STY_;
            RME(i) = RME_;
            
            f = [-log(STY) -log(RME)]; % maximise STY, maximise RME
            
% Update data set
                
            data.x(i,:) = conditions(i,:);
            data.y(i,:) = f(i,:);  
            
% Determine number of Pareto optimal solutions after each iteration
            
            STY_actual = exp(-data.y(:,1));
            RME_actual = exp(-data.y(:,2));

            data.z = [STY_actual RME_actual];
            [opt, data.idxs] = paretoFront(data.z);
            data.optconds = data.x(data.idxs,:); % extract conditions corresponding to Pareto optimal solutions
            opt_size = size(opt);
            opt_n = opt_size(1,1);
            data.opt(i,:) = opt_n;
            data.front = opt; 
            
% Calculate hypervolume after each iteration (requires known reference points)

            F = opt; 
            AU = [566.7668 3.0476];  
            U = [1.1149e+03 6.8251];
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
            
            tR = conditions(:,1); % in min 
            T = 50; % in oC, rate function converts
            ConcA = 1;
            equivB = conditions(:,2);
            ConcB = ConcA.*equivB;
            
            SM_yield = zeros((size(conditions(:,1),1)) , 1);
            SM2_yield = zeros((size(conditions(:,1),1)) , 1);
            intermediate_yield = zeros((size(conditions(:,1),1)) , 1);
            product_yield = zeros((size(conditions(:,1),1)) , 1);
            E_factor = zeros((size(conditions(:,1),1)) , 1);
            STY = zeros((size(conditions(:,1),1)) , 1);
            RME = zeros((size(conditions(:,1),1)) , 1);
            
for i=1:1:(size(conditions(:,1),1))
    
            [Cfinal, Yield, TimeData, ConcData] = PK_Constants(ConcA, ConcB(i), T * ones(1,4), tR(i));
            
            % Calculate percentage yields
            
            SM = Yield(1)*100;
            SM2 = Yield(2)*100;
            intermediate = Yield(3)*100;
            product = Yield(4)*100;
            
            % Adjust yields with error
            
            SM = SM+((rand-0.5)*0.5)+(SM*((rand-0.5)/100));
            if SM < 0
                SM = 1*10^-6;
            end
            if SM > 100
                SM = 100;
            end
            
            SM2 = SM2+((rand-0.5)*0.5)+(SM2*((rand-0.5)/100));
            if SM2 < 0
                SM2 = 1*10^-6;
            end
            if SM2 > 100
                SM2 = 100;
            end
            
            intermediate = intermediate+((rand-0.5)*0.5)+(intermediate*((rand-0.5)/100));
            if intermediate < 0
                intermediate = 1*10^-6;
            end
            if intermediate > 100
                intermediate = 100;
            end 
            
            product = product+((rand-0.5)*0.5)+(product*((rand-0.5)/100));
            if product < 0
                product = 1*10^-6;
            end
            if product > 100
                product = 100;
            end
            
            % Calculate other objectives
            
            STY_ = STY_calc(product/100, tR(i), ConcData(1,1), 139.20);

            RME_ = (139.20*product)/(114.14+(61.08*(ConcData(1,5+1)./ConcData(1,1)))); % 5 = global 'n' (see Lactose constants)
            
            % Update data
            
            SM_yield(i) = SM;
            SM2_yield(i) = SM2;
            intermediate_yield(i) = intermediate;
            product_yield(i) = product;
            STY(i) = STY_;
            RME(i) = RME_;

            f = [-log(STY) -log(RME)]; % maximise STY, maximise RME
            
% Update data set
                
            data.x = [data.x ; conditions(i,:)];
            data.y = [data.y ; f(i,:)];
            total_expts = (size(data.x,1)); % calculate total experiments so far       
            
% Determine number of Pareto optimal solutions after each iteration
            
            STY_actual = exp(-data.y(:,1));
            RME_actual = exp(-data.y(:,2));

            data.z = [STY_actual RME_actual];
            [opt, data.idxs] = paretoFront(data.z);
            data.optconds = data.x(data.idxs,:); % extract conditions corresponding to Pareto optimal solutions
            opt_size = size(opt);
            opt_n = opt_size(1,1);
            data.opt = [data.opt ; opt_n];
            data.front = opt;       
            
% Calculate hypervolume after each iteration (requires known reference points)

            F = opt; 
            AU = [566.7668 3.0476];  
            U = [1.1149e+03 6.8251];
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