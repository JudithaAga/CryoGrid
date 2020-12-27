theta_o= 0;

for theta_m = [0.05:0.05:0.95]
    
    porosity = 1 - theta_o - theta_m;
    
    c_w = 4.2e6;
    c_i = 1.9e6;
    c_m = 2e6;
    c_o = 2.5e6;
    L_sl = 3.34e8;
    T0=273.15;
    beta_interface = 1./3;
    
    T_min = -60;
    sat_waterIce_min = 0.005;
    
    alpha = 8e-4; %(Pa^-1) this is not the same alpha as we have used, it is the conversion form 1/m to 1Pa (factor 10-4) 
    m=0.19;
    
    
    alpha = 1.11e-4;
    n=1.48;
%     alpha = 4e-4;
%     n=2;
%     alpha = 1.49e-4;
%     n=1.25;
    m=1-1./n;
    
    % sat_waterIce_list = [0.2; 0.5; 1];
    T=[-60:0.01:0]';
    
      n=1./(1-m);  

    
    % for i=1:3
    %     sat_waterIce = sat_waterIce_list(i)
    sat_waterIce_list = [0.005;0.05; 0.5; 0.8; 1] %[0.01:0.01:1];
    for i=1:size(sat_waterIce_list,1)  
        
        sat_waterIce = sat_waterIce_list(i);
        %sat_waterIce=0.05;
        
        %matrix water pressure in unfrozen state
        mwp0 = 1./alpha .* ((sat_waterIce.^(-1./m)-1)).^(1./n);
        
        %matrix water potential
        mpw = -L_sl.*T./T0 .*beta_interface .* double(T<0) + mwp0;
        
        %freeze curve
        sat_water = double(mpw>0) .* (1+(alpha.*mpw).^n).^(-m) + double(mpw<=0);
        sat_ice = sat_waterIce-sat_water;
        % plot(T,sat_water)
        % hold on
        
        %energy
        energy = T.* (theta_m .* c_m + theta_o .* c_o + porosity .* sat_waterIce .* (c_w .* double(T>=0)+ c_i.* double(T<0)));
        energy = energy - double(T<0) .* L_sl .* porosity .* (sat_waterIce - sat_water);
        
        %   plot(energy, T)
        
        
        C0 = theta_m .* c_m + theta_o .* c_o + porosity .* sat_waterIce .* c_i;
        X = -L_sl./ T0 .* beta_interface;
        L0 = porosity.* L_sl;
        
        %dimensionless quantities computed from the "true" quantities
        E_prime = energy ./ L0 + sat_waterIce + C0 .* mwp0 ./ X./ L0; %conversion of energy to E_prime (dimensionless)
        T_prime = alpha .* (X .* T + mwp0); %conversion of temperature to T_prime (dimensionless)
        
        
        %This is the second parameter in dimensionless equation 
        gamma = C0 ./ alpha ./ X ./ L0 ;
        
        T_acc = 0;
        
%         for ind=-1:2:1
%             
%             gamma = gamma0 ;%+ ind.*3.7e-6 /100 ;
        %dimensionless quantities computed from the dimensionless equation

        T_prime_max = alpha .* X .* T_min + ((sat_waterIce_min.^(-1./m)-1)).^(1./n);
        E_prime_min = gamma .* T_prime_max + (1 + T_prime_max.^n).^(-m);
        
        T_prime_max = T_prime_max ;
        E_prime_min = E_prime_min ;
        
        %start with 0 and T_prime_max divided in something like 2^10 segments,
        %then subdivide them by two, calculate the interpolated and the real
        %point, and throw out the 2^10 points for which the interpolation
        %performs best
        %, 
%         LUT_size = 2.^10;
%         
%         target_delta_E = -1;%-(1 - E_prime_min) ./ LUT_size;
%         
% 
%         T_prime2 = [1e-1; 0];
%         E_prime2 = [gamma .* T_prime2(1) + (1 + T_prime2(1).^n).^(-m); 1];
%          while E_prime2(1) > E_prime_min
%            dE_dT = gamma - m.*(1 + T_prime(1).^n).^(-m-1) .* n .* T_prime(1).^(n-1);
%            T_prime2_projected = T_prime2(1) + target_delta_E ./ dE_dT;
%            E_prime2_projected = gamma .* T_prime2_projected + (1 + T_prime2_projected.^n).^(-m);
%            %residual = (E_prime2_projected - E_prime2(1)) ./ target_delta_E;
%            deltaE_deltaT = (E_prime2_projected - E_prime2(1)) ./ (T_prime2_projected - T_prime2(1));
%            tolerance = abs((dE_dT - deltaE_deltaT) ./ dE_dT);
%            if tolerance > 0.01
%                target_delta_E = target_delta_E ./2;
%            elseif tolerance < 0.001
%                target_delta_E = target_delta_E .*2;
%            else
%             T_prime2 = [T_prime2_projected; T_prime2];
%             E_prime2 = [E_prime2_projected; E_prime2];
%             wjdfjhwe
%            end
%             
%          end
        
        LUT_size = 2.^10;

        T_prime2 = linspace(T_prime_max, 0, LUT_size+1)';
        E_prime2 = gamma .* T_prime2 + (1 + T_prime2.^n).^(-m);


        for count = 1:15
%             plot(E_prime2)
%             hold on
            %add points
            add_points_T =[(T_prime2(1:end-1,1) + T_prime2(2:end,1))./2  (3.*T_prime2(1:end-1,1) + T_prime2(2:end,1))./4 (T_prime2(1:end-1,1) + 3.*T_prime2(2:end,1))./4];%midpoints in T
            add_points_E = gamma .* add_points_T + (1 + add_points_T.^n).^(-m);
            
            delta_E_prime2 = repmat(E_prime2(2:end,1) - E_prime2(1:end-1,1), 1, 3) ;
            delta_T_prime2 = repmat(T_prime2(2:end,1) - T_prime2(1:end-1,1), 1, 3);
            add_points_T_interp = repmat(T_prime2(1:end-1),1,3) + (add_points_E - repmat(E_prime2(1:end-1),1,3)) ./ delta_E_prime2 .* delta_T_prime2;
            
            add_points_T = add_points_T(:);
            add_points_E = add_points_E(:);
            mismatch = abs(add_points_T_interp(:) - add_points_T);
            [~,sort_index]=sort(mismatch);
            sort_index=sort_index(end-size(sort_index,1)./3+1:end,1);
            T_prime2 = [T_prime2; add_points_T(sort_index,1)];
            E_prime2 = [E_prime2; add_points_E(sort_index,1)];
            [T_prime2, sort_index] = sort(T_prime2,'descend');
            E_prime2=E_prime2(sort_index,1);
            
            %reduce points
            reduce_points_T = [T_prime2(2:4:end-3,1) T_prime2(3:4:end-2,1) T_prime2(4:4:end-1,1)];
            reduce_points_E = [E_prime2(2:4:end-3,1) E_prime2(3:4:end-2,1) E_prime2(4:4:end-1,1)];
            T_prime2 = T_prime2(1:4:end,1);
            E_prime2 = E_prime2(1:4:end,1);
            
            delta_E_prime2 = repmat(E_prime2(2:end,1) - E_prime2(1:end-1,1), 1, 3) ;
            delta_T_prime2 = repmat(T_prime2(2:end,1) - T_prime2(1:end-1,1), 1, 3);
            reduce_points_T_interp = repmat(T_prime2(1:end-1),1,3) + (reduce_points_E - repmat(E_prime2(1:end-1),1,3)) ./ delta_E_prime2 .* delta_T_prime2;
            
            reduce_points_T = reduce_points_T(:);
            reduce_points_E = reduce_points_E(:);
            mismatch = abs(reduce_points_T_interp(:) - reduce_points_T);
            [~,sort_index]=sort(mismatch);
            sort_index=sort_index(end-size(sort_index,1)./3+1:end,1);
            T_prime2 = [T_prime2; reduce_points_T(sort_index,1)];
            E_prime2 = [E_prime2; reduce_points_E(sort_index,1)];
            [T_prime2, sort_index] = sort(T_prime2,'descend');
            E_prime2 = E_prime2(sort_index,1);
            
        end

        
%         plot(E_prime2, T_prime2,'+')

% 
%         
%                 real_points = T_prime2(2:2:end-1,1);
%                 mismatch = abs(interp_points - real_points);
%                 [~,sort_index]=sort(mismatch);
%                 max_error=[max_error; max(mismatch)];
%                 %reduce
%                 %delete_index = sort_index(end-LUT_size./4+1:end,1) .* 2;
%                 delete_index = sort_index(1:LUT_size./4,1) .* 2;
%                 T_prime2(delete_index,:) = [];
%                 E_prime2 = gamma .* T_prime2 + (1 + T_prime2.^n).^(-m);
%             %end
%             %expand
%             T_prime2_fill = (T_prime2(1:end-1,1)+T_prime2(2:end,1))./2;
%             E_prime2_fill = gamma .* T_prime2_fill + (1 + T_prime2_fill.^n).^(-m);
%             E_prime2_new = zeros(LUT_size+1,1);
%             T_prime2_new = zeros(LUT_size+1,1);
%             
%             E_prime2_new(1:2:end,1) = E_prime2;
%             E_prime2_new(2:2:end-1,1) = E_prime2_fill;
%             T_prime2_new(1:2:end,1) = T_prime2;
%             T_prime2_new(2:2:end-1,1) = T_prime2_fill;
%             
%             T_prime2 = T_prime2_new;
%             E_prime2 = E_prime2_new;

%             plot(E_prime2)
%             hold on
        %end
%         max_error=[];
%         
%         for ind=1:20
%             for iter = 1:2
%                 interp_points_T =(T_prime2(1:2:end-2,1) + T_prime2(3:2:end,1))./2;
%                 real_points = T_prime2(2:2:end-1,1);
%                 mismatch = abs(interp_points - real_points);
%                 [~,sort_index]=sort(mismatch);
%                 max_error=[max_error; max(mismatch)];
%                 %reduce
%                 %delete_index = sort_index(end-LUT_size./4+1:end,1) .* 2;
%                 delete_index = sort_index(1:LUT_size./4,1) .* 2;
%                 T_prime2(delete_index,:) = [];
%                 E_prime2 = gamma .* T_prime2 + (1 + T_prime2.^n).^(-m);
%             end
%             %expand
%             T_prime2_fill = (T_prime2(1:end-1,1)+T_prime2(2:end,1))./2;
%             E_prime2_fill = gamma .* T_prime2_fill + (1 + T_prime2_fill.^n).^(-m);
%             E_prime2_new = zeros(LUT_size+1,1);
%             T_prime2_new = zeros(LUT_size+1,1);
%             
%             E_prime2_new(1:2:end,1) = E_prime2;
%             E_prime2_new(2:2:end-1,1) = E_prime2_fill;
%             T_prime2_new(1:2:end,1) = T_prime2;
%             T_prime2_new(2:2:end-1,1) = T_prime2_fill;
%             
%             T_prime2 = T_prime2_new;
%             E_prime2 = E_prime2_new;
% 
% %             plot(E_prime2)
% %             hold on
%         end
        
        
%         LUT_size = 2.^18+1;
%         E_prime3 = linspace(E_prime2_min, 1, LUT_size);
%         
%         T_prime3 = interp1(E_prime2, T_prime2, E_prime3);


        
%         curvature = -m.*(-m-1).*(1+T_prime2.^n).^(-m-2) .* n.* T_prime2.^(n-1) .* n.* T_prime2.^(n-1) + -m.*(1+T_prime2.^n).^(-m-1).*n.*(n-1).*T_prime2.^(n-2);
%         index=size(curvature,1);
%         while abs(curvature(index,1))<2.5e-12
%             index=index-1;
%         end
%         index2=size(curvature,1);
%         while abs(curvature(index2,1))<2.5e-9
%             index2=index2-1;
%         end
% 
%         E_range1 = E_prime2(index);
%         E_range2 = E_prime2(index2);
%         LUT_size=2.^16;
%         if E_range1~=E_prime2_min
%             E_prime3 = linspace(E_prime2_min, E_range1, LUT_size+1);
%             E_prime3 = E_prime3(1:end-1);
%         else
%             E_prime3 = [];
%         end
%         E_prime3 = [E_prime3 linspace(E_range1, E_range2, LUT_size)];
%         E_prime3 = E_prime3(1:end-1);        
%         E_prime3 = [E_prime3 linspace(E_range2, 1, 2.*LUT_size)];
%         T_prime3 = interp1(E_prime2, T_prime2, E_prime3);
        
        
        
        %Plot the two on top of each other in "prime space"
%         plot(E_prime2, T_prime2)
%         
%         % figure
% %          on
%         plot(E_prime3, T_prime3)
%         
%         plot(E_prime2, T_prime2, '+')
%         hold on
%         plot(E_prime, T_prime)
%         

        %---------------interpolate and test
        test_T = [-50:0.01:-0.01];
        
        mpw_test = -L_sl.*test_T./T0 .*beta_interface .* double(test_T<0) + mwp0;
        sat_water_test = double(mpw_test>0) .* (1+(alpha.*mpw_test).^n).^(-m) + double(mpw_test<=0);
        sat_ice_test = sat_waterIce - sat_water_test;
        energy_test = test_T.* (theta_m .* c_m + theta_o .* c_o + porosity .* sat_waterIce .* (c_w .* double(test_T>=0)+ c_i.* double(test_T<0)));
        energy_test = energy_test - double(test_T<0) .* L_sl .* porosity .* (sat_waterIce - sat_water_test);
        E_prime_test = energy_test ./ L0 + sat_waterIce + C0 .* mwp0 ./ X./ L0; 
        
%         pos = (E_prime_test - 1) ./  (E_prime3(1)- E_prime3(2));
        T_prime_test = interp1(E_prime2, T_prime2, E_prime_test, 'linear');
        
        T_test = (T_prime_test ./ alpha - mwp0) ./ X;

%         T_acc = T_acc + T_test;
%         end
        
         plot(test_T, test_T - T_test)
         hold on
        
    end
end

% porosity_max = 0.95;
% gamma_max =  ((1-porosity_max).* c_o ) ./ alpha ./ X ./ (porosity_max.* L_sl);
% log(-gamma_max)
% porosity_min = 0.05;
% gamma_min =  ((1-porosity_min).* c_m  + c_i) ./ alpha ./ X ./ (porosity_min.* L_sl);
% log(-gamma_min)


%




% end

%  plot(energy, sat_water)



% B = sat_water./(1-sat_ice);
% P_cgl = 1./alpha .* ((B.^(-1./m)-1)).^(1./n);
% P_csl = -L_sl.*T./T0 .*double(T<0);

% get_E(-10)
%
%
%
% fsolve(@get_E, -15)
% 
% function energy = get_E(T)
% 
% porosity = 0.5;
% alpha = 8e-4; %(Pa^-1)
% m=0.19;
% L_sl = 3.34e8;
% T0=273.15;
% 
% alpha = 1.11e-4;
% n=1.48;
% m=1-1./n;
% 
% % sat_waterIce_list = [0.2; 0.5; 1];
% 
% 
% n=1./(1-m);
% 
% % for i=1:3
% %     sat_waterIce = sat_waterIce_list(i)
% 
% sat_waterIce = 0.5;
% A = -L_sl.*T./T0 ./3 .*double(T<0) + 1./alpha .* ((sat_waterIce.^(-1./m)-1)).^(1./n);
% 
% sat_water = double(A>0) .* (1+(alpha.*A).^n).^(-m) + double(A<=0);
% sat_ice = sat_waterIce-sat_water;
% 
% 
% energy = T.* (porosity .* 2e6 + (1-porosity) .* sat_waterIce .* (4.2e6 .* double(T>=0)+ 1.9e6.* double(T<0)));
% energy = energy - double(T<0) .* 3.34e8 .* (1-porosity) .* sat_ice;
% 
% energy= energy +8.9414e+07;
% end