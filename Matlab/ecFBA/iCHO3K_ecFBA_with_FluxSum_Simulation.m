load('Context_Specific_Models.mat');
Models = struct2cell(Context_Specific_Models);

pick = input(['Choose the Data\n' ...
                '   [1] WT\n' ...
                '   [2] ZeLa\n' ...
                'Response: ']);

switch pick
    case 1
        Model_id = {'WT_P4_U1';'WT_P4_U2';'WT_P4_U3';'WT_P6_U1';'WT_P6_U2';'WT_P6_U3'};
        WT_rates = readtable('uptake_secretion_rates_wt_sorted.csv');
        Rxn_id = table2array(WT_rates(:,1));
        Rates = WT_rates(:,[4,10,16,5,11,17]);
        Rates.Properties.RowNames = Rxn_id;
        Rates.Properties.VariableNames = Model_id;
        grRate = table2array(Rates('biomass_cho_s',:));
        Rates('biomass_cho_s',:) = [];
        Rxn_id = Rates.Properties.RowNames;
        Models = Models(1:6,1);

    case 2
        Model_id = {'ZeLa_P4_U4';'ZeLa_P4_U5';'ZeLa_P4_U6';'ZeLa_P4_U7';'ZeLa_P4_U8';...
                    'ZeLa_P6_U4';'ZeLa_P6_U6';'ZeLa_P6_U7';'ZeLa_P6_U8'};
        ZeLa_rates = readtable('uptake_secretion_rates_zela_sorted.csv');
        Rxn_id = table2array(ZeLa_rates(:,1));
        Rates = ZeLa_rates(:,[4,10,16,22,28,5,17,23,29]);
        Rates.Properties.RowNames = Rxn_id;
        Rates.Properties.VariableNames = Model_id;
        grRate = table2array(Rates('biomass_cho_s',:));
        Rates('biomass_cho_s',:) = [];
        Rxn_id = Rates.Properties.RowNames;
        Models = Models(7:end,1);
    otherwise
        disp('Invalid Input :(');
end

sample_size = input('Enter Sample Size Required: ');

modelIrrev = struct();
Solutions = struct();
fluxSum = struct();
        for i = 1:length(Models)
            model_rate = table2array(Rates(:,i));
            test_model = Models{i,1};
            fprintf('\n%s loaded\n',string(Model_id(i)));
            rev = zeros(length(test_model.rxns),1);
            for r = 1:length(test_model.rxns)
                if test_model.lb(r)<0 && test_model.ub(r)>0
                    rev(r,1) = 1;
                end
            end
            test_model.rev = rev;
            disp('Reversablility Updated');
            test_model.kcat_f = str2double(test_model.kcat_f);
            test_model.kcat_f(find(isnan(test_model.kcat_f))) = 0;
            test_model.kcat_b = str2double(test_model.kcat_b);
            test_model.kcat_b(find(isnan(test_model.kcat_b))) = 0;
            test_model.molwt = str2double(test_model.molwt);
            test_model.molwt(find(isnan(test_model.molwt))) = 0;

            test_model = changeRxnBounds(test_model,{'SK_pre_prot_r','SK_Ser_Thr_g','SK_Tyr_ggn_c','SK_Asn_X_Ser_Thr_r'},-0.1,'l');
            test_model = changeRxnBounds(test_model,{'SK_pre_prot_r','SK_Ser_Thr_g','SK_Tyr_ggn_c','SK_Asn_X_Ser_Thr_r'},1000,'u');
            Unbound_rxns = {'EX_ca2_e';'EX_cl_e';'EX_co_e';'EX_co2_e';'EX_fe2_e';'EX_h_e';'EX_h2o_e';...
                            'EX_hco3_e';'EX_i_e';'EX_k_e';'EX_na1_e';'EX_nh4_e';'EX_no_e';'EX_o2_e';...
                            'EX_pi_e';'EX_so4_e';'EX_etoh_e'};
            test_model = changeRxnBounds(test_model,Unbound_rxns,-1000,'l');
            test_model = changeRxnBounds(test_model,Unbound_rxns,1000,'u');
            disp('Inorganics Exchange Reactions Relaxed');
            test_model = changeRxnBounds(test_model,'GLDBRAN',0,'b');
            test_model = changeRxnBounds(test_model,'EHGLAT2m',0,'b');
            disp('Loop Reactions GLDBRAN and EHGLAT2m set to Zero');            
            disp('Uptake AA Rates Constrained with 10% error on lower side & Secretory AA Rates constrained on upper side...');
            for j = 1:length(Rxn_id)
                if model_rate(j)<0
                    test_model = changeRxnBounds(test_model,Rxn_id(j),model_rate(j)*1.1,'l');
                    test_model = changeRxnBounds(test_model,Rxn_id(j),0,'u');
                else
                    test_model = changeRxnBounds(test_model,Rxn_id(j),0,'l');
                    test_model = changeRxnBounds(test_model,Rxn_id(j),model_rate(j)*1.1,'u');
                end
            end

                switch pick
                    case 1
                        test_model = changeRxnBounds(test_model,'EX_glc_e',model_rate(13)*1.1,'l');
                        test_model = changeRxnBounds(test_model,'EX_glc_e',model_rate(13)*0.9,'u');
                        disp("Glucose constrained with 10% error");
                        test_model = changeRxnBounds(test_model,'biomass_cho_s',grRate(i)*0.9,'l');
                        test_model = changeRxnBounds(test_model,'biomass_cho_s',grRate(i)*1.1,'u');
                        disp("Biomass constrained with 10% error");

                        test_model = changeRxnBounds(test_model,'EX_lac_L_e',model_rate(20)*0.9,'l');
                        test_model = changeRxnBounds(test_model,'EX_lac_L_e',model_rate(20)*1.1,'u');
                        disp("l-Lactate uptake unconstrained and secretion constrained w 10% error");
                    case 2
                        test_model = changeRxnBounds(test_model,'EX_glc_e',model_rate(13)*1.1,'l');
                        test_model = changeRxnBounds(test_model,'EX_glc_e',model_rate(13)*0.9,'u');
                        disp("Glucose constrained with 10% error");
                        test_model = changeRxnBounds(test_model,'biomass_cho_s',grRate(i)*0.9,'l');
                        test_model = changeRxnBounds(test_model,'biomass_cho_s',grRate(i)*1.1,'u');
                        disp("Biomass constrained with 10% error");
                        test_model = changeRxnBounds(test_model,'LDH_L',0,'b');
                        disp("LDH_L Knocked out");

                        if model_rate(21) ~=0
                            test_model = changeRxnBounds(test_model,'EX_lac_L_e',model_rate(20)*0.9,'l');
                            test_model = changeRxnBounds(test_model,'EX_lac_L_e',model_rate(20)*1.1,'u');
                            disp("l-Lactate uptake and secretion constrained w 10% error");
                        else
                            test_model = changeRxnBounds(test_model,'EX_lac_L_e',0,'b');
                            disp("l-Lactate hard bound to 0");
                        end
                    otherwise
                        disp('Invalid Input :(');
                end
                fprintf('\n%s Minimization ecFBA Simulation Initiated ...\n',string(Model_id(i)));
                [modelIrrev.(strcat("modelIrrev_",Model_id(i))), Solutions.(strcat("ec_Solutions_",Model_id(i))),fluxSum.(strcat("fluxSum_",Model_id(i)))] = minEcFBAwithFluxSum(test_model,'biomass_cho_s',sample_size);
                
        end
