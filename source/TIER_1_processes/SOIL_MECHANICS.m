%========================================================================
% CryoGrid TIER1 library class, functions related to soil mechanics
% J. Schmidt, May 2021
%========================================================================

classdef SOIL_MECHANICS < BASE

    methods
        
        %Calculate settlement out of water flow, as well as overburden pressure and bearing capacity for the next timestep
        function ground = soil_mechanics(ground, tile)
            
            %Porosity, void ratio and saturation is updated as this changed for saturated cells in advanced prognostics
            ground.STATVAR.Xporosity = 1 - (ground.STATVAR.mineral + ground.STATVAR.organic) ./ ground.STATVAR.layerThick ./ ground.STATVAR.area;
            ground.STATVAR.porosity = min(ground.STATVAR.Xporosity, ground.STATVAR.initial_porosity);
            ground.STATVAR.Xvoid_ratio = ground.STATVAR.Xporosity ./ (1 - ground.STATVAR.Xporosity);
            ground.STATVAR.void_ratio = ground.STATVAR.porosity ./ (1 - ground.STATVAR.porosity);
            ground.STATVAR.saturation = (ground.STATVAR.waterIce ./ (ground.STATVAR.layerThick .* ground.STATVAR.area)) ./ ground.STATVAR.Xporosity;
            
            %Find all unsaturated cells --> will be updated in this function 
            unsaturated = ground.STATVAR.saturation <= 1-1e-6;

            %Calculate the overburden pressure for each cell
            %Below saturation threshold (cell is saturated less than 50%), all overburden pressure is on the soil
            %Above saturation threshold (cell is saturated 50% or more), overburden pressure is reduced linearly with increasing saturation
            %At 100% saturation, overburden is reduced by density of water mineral content * (density mineral - density water)
%             threshold = 0.5; above_threshold = ground.STATVAR.saturation >= threshold;
%             overburden_pressure_per_cell_normal = (ground.STATVAR.mineral .* ground.CONST.rho_m + ground.STATVAR.organic .* ground.CONST.rho_o + ...
%                 ground.STATVAR.ice .* ground.CONST.rho_i + ground.STATVAR.water .* ground.CONST.rho_w) ./ ground.STATVAR.area .*ground.CONST.g; %[Pa]
%             overburden_pressure_per_cell_normal_boyant = (ground.STATVAR.mineral .* (ground.CONST.rho_m - ground.CONST.rho_w) + ground.STATVAR.organic .* (ground.CONST.rho_o - ground.CONST.rho_w) + ...
%                 ground.STATVAR.ice .* (ground.CONST.rho_i - ground.CONST.rho_w) + (ground.STATVAR.layerThick - ground.STATVAR.mineral - ground.STATVAR.organic - ground.STATVAR.ice - ground.STATVAR.water) .* (-ground.CONST.rho_w)) ./ ground.STATVAR.area .*ground.CONST.g; %[Pa]
%             overburden_pressure_per_cell = double(~above_threshold) .* overburden_pressure_per_cell_normal + ...
%                 double(above_threshold) .* (overburden_pressure_per_cell_normal + (ground.STATVAR.saturation-threshold)./ (1-threshold) .* (overburden_pressure_per_cell_normal_boyant - overburden_pressure_per_cell_normal));
            threshold = 0.5; above_threshold = ground.STATVAR.saturation >= threshold;
            overburden_pressure_per_cell_normal = (ground.STATVAR.mineral .* ground.CONST.rho_m + ground.STATVAR.organic .* ground.CONST.rho_o + ...
                ground.STATVAR.waterIce .* ground.CONST.rho_w) ./ ground.STATVAR.area .*ground.CONST.g; %[Pa]
            overburden_pressure_per_cell_normal_boyant = (ground.STATVAR.mineral .* (ground.CONST.rho_m - ground.CONST.rho_w) + ground.STATVAR.organic .* (ground.CONST.rho_o - ground.CONST.rho_w) + ...
                + (ground.STATVAR.layerThick - ground.STATVAR.mineral - ground.STATVAR.organic - ground.STATVAR.waterIce) .* (-ground.CONST.rho_w)) ./ ground.STATVAR.area .*ground.CONST.g; %[Pa]
            overburden_pressure_per_cell = double(~above_threshold) .* overburden_pressure_per_cell_normal + ...
                double(above_threshold) .* (overburden_pressure_per_cell_normal + (ground.STATVAR.saturation-threshold)./ (1-threshold) .* (overburden_pressure_per_cell_normal_boyant - overburden_pressure_per_cell_normal));
  
            obp_old = ground.STATVAR.overburden_pressure;
            ground.STATVAR.overburden_pressure = cumsum(overburden_pressure_per_cell)-overburden_pressure_per_cell./2; %[Pa]
            
            %Add external pressure
            ground.STATVAR.overburden_pressure = ground.STATVAR.overburden_pressure + ground.PARA.external_pressure;
            
            %Make mean of current overburden pressure and n times old overburden pressure --> needed for stable calculations
            ground.STATVAR.overburden_pressure = (ground.STATVAR.overburden_pressure + ((ground.PARA.smoothing_factor - 1).*obp_old))./ground.PARA.smoothing_factor; 
            
            %Update porosity and layerThick for unsaturated cells: not Xporosity because unsaturated gridcells can't have a higher porosity than initial_porosity
            ground.STATVAR.porosity(unsaturated) = (ground.STATVAR.initial_voidRatio(unsaturated) - ground.STATVAR.compression_index(unsaturated) .* log10(ground.STATVAR.overburden_pressure(unsaturated)./ground.STATVAR.reference_pressure(unsaturated))) ./ ...
                (1 + ground.STATVAR.initial_voidRatio(unsaturated) - ground.STATVAR.compression_index(unsaturated) .* log10(ground.STATVAR.overburden_pressure(unsaturated)./ground.STATVAR.reference_pressure(unsaturated)));
            %do not allow higher porosities than initial void ration
            ground.STATVAR.porosity(unsaturated) = min(ground.STATVAR.porosity(unsaturated), ground.STATVAR.initial_porosity(unsaturated));
            ground.STATVAR.layerThick(unsaturated) = max(((ground.STATVAR.mineral(unsaturated) + ground.STATVAR.organic(unsaturated)) ./ (1 - ground.STATVAR.porosity(unsaturated))) ./ ground.STATVAR.area(unsaturated), ...
                (ground.STATVAR.mineral(unsaturated) + ground.STATVAR.organic(unsaturated)+ ground.STATVAR.waterIce(unsaturated))./ ground.STATVAR.area(unsaturated));
  
            %Recalculate porosity, void ratio and saturation for all cells            
            ground.STATVAR.Xporosity = 1 - (ground.STATVAR.mineral + ground.STATVAR.organic) ./ ground.STATVAR.layerThick ./ ground.STATVAR.area;
            ground.STATVAR.porosity = min(ground.STATVAR.Xporosity, ground.STATVAR.initial_porosity);
            ground.STATVAR.saturation = (ground.STATVAR.waterIce ./ (ground.STATVAR.layerThick .* ground.STATVAR.area)) ./ ground.STATVAR.Xporosity;
%            ground.STATVAR.water_saturation = ((ground.STATVAR.water + ground.STATVAR.Xwater) ./ (ground.STATVAR.layerThick .* ground.STATVAR.area)) ./ ground.STATVAR.porosity;
            ground.STATVAR.Xvoid_ratio = ground.STATVAR.Xporosity ./ (1 - ground.STATVAR.Xporosity);
            ground.STATVAR.void_ratio = ground.STATVAR.porosity ./ (1 - ground.STATVAR.porosity);
            ground.STATVAR.layerThick = ((ground.STATVAR.mineral + ground.STATVAR.organic) ./ (1 - ground.STATVAR.Xporosity)) ./ ground.STATVAR.area;

            %Calculate bearing capacity (= effective stress on soil)
            bc_old = ground.STATVAR.bearing_capacity;
            ground.STATVAR.bearing_capacity = (10.^((ground.STATVAR.initial_voidRatio - ground.STATVAR.Xporosity - ground.STATVAR.initial_voidRatio .* ground.STATVAR.Xporosity) ./ (ground.STATVAR.compression_index - ground.STATVAR.Xporosity .* ground.STATVAR.compression_index))) .* ground.STATVAR.reference_pressure;
            ground.STATVAR.bearing_capacity = max(ground.STATVAR.reference_pressure, ground.STATVAR.bearing_capacity); %generate "Xice" (void ratio becomes smaller than initial void ratio)
            ground.STATVAR.bearing_capacity = (ground.STATVAR.bearing_capacity + (ground.PARA.smoothing_factor - 1) .* bc_old) ./ ground.PARA.smoothing_factor;  
            
            %Update waterIce
            waterIce = ground.STATVAR.saturation .* ground.STATVAR.Xporosity;
            ground.STATVAR.waterIce = waterIce .* ground.STATVAR.layerThick .* ground.STATVAR.area;
            XwaterIce = ground.STATVAR.saturation .* (ground.STATVAR.Xporosity - ground.STATVAR.porosity);
            ground.STATVAR.XwaterIce = XwaterIce .* ground.STATVAR.layerThick .* ground.STATVAR.area;
%            ground.STATVAR.XwaterIce = max(0,ground.STATVAR.waterIce - ground.STATVAR.saturation .* ground.STATVAR.porosity .* ground.STATVAR.initial_layerThick .* ground.STATVAR.area);
%    
            
            %Split waterIce up in XwaterIce and waterIce
%            porosity_equilibrium = ground.STATVAR.initial_voidRatio ./ (1 + ground.STATVAR.initial_voidRatio);
%            waterIce_equilibrium = (ground.STATVAR.mineral + ground.STATVAR.organic) ./ (1 - porosity_equilibrium) - ground.STATVAR.mineral - ground.STATVAR.organic;
            
%            ground.STATVAR.XwaterIce = ground.STATVAR.waterIce.*0;
%            excessIce = ground.STATVAR.saturation >=1-1e-6 & waterIce_equilibrium < ground.STATVAR.waterIce;
%            ground.STATVAR.XwaterIce(excessIce) = ground.STATVAR.waterIce(excessIce) - waterIce_equilibrium(excessIce);
%            ground.STATVAR.waterIce(excessIce) = waterIce_equilibrium(excessIce);
        end
    end
end

