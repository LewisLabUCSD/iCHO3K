function [dFBA_results] = COSMIC_dFBA(model,time_range,dFBA_data,fix_flx,fix_comp,notes)
wbar = waitbar(0,'Initializing dFBA Simulation');
np = size(dFBA_data.components,1);
[time,phase_transition] = deal(zeros(10000,1));
[profiles,flux_growth,flux_prod] = deal(zeros(10000,np));
ctr = 1;
rtol = 1e-2;
h = 1e-4;
m1 = model;
t_curr = time_range(1);
done = false;
c0 = dFBA_data.initial_concentrations(:);
c0 = c0';
cmp_names = dFBA_data.component_names;
components = dFBA_data.components;
rxns = dFBA_data.components;
time(ctr) = t_curr;
profiles(ctr,:) = c0;
c_in = zeros(size(profiles,2),1);
for i = 1:length(c_in)
    c_in(i) = dFBA_data.Perfusion.concentrations(ismember(dFBA_data.Perfusion.components,cmp_names(i)));
end
c_in = c_in';
kinetic_comps = dFBA_data.kinetic.components;
kinetic_growth = [dFBA_data.kinetic.Vm_growth(:),dFBA_data.kinetic.Km(:)]';
kinetic_prod = [dFBA_data.kinetic.Vm_prod(:),dFBA_data.kinetic.Km(:)]';
p_growth = dFBA_data.Objectives.priority_growth;
f_growth = dFBA_data.Objectives.c_growth;
msk = f_growth>=0;
p_growth = p_growth(msk);
f_growth = f_growth(msk);
p_prod = dFBA_data.Objectives.priority_prod;
f_prod = dFBA_data.Objectives.c_prod;
msk = f_prod >= 0;
p_prod = p_prod(msk);
f_prod = f_prod(msk);
classification_data = dFBA_data.Phase_classification;

while ~done
    
    % Single Step
    kc = zeros(size(kinetic_comps));
    for i = 1:length(kc)
        kc(i) = c0(ismember(components,kinetic_comps(i)));
    end
    vkg = evaluate_fluxes(kc,kinetic_growth);
    vkp = evaluate_fluxes(kc,kinetic_prod);
    vg = simulate_model(model,components,kinetic_comps,vkg,p_growth,f_growth,fix_flx,fix_comp);
    vp = simulate_model(model,components,kinetic_comps,vkp,p_prod,f_prod,fix_flx,fix_comp);
    
    f = phase_progress(c0,components,classification_data);
    v = (1-f)*vg + f*vp;
    dc = eval_reactor(c0,v,c_in);
    c1 = c0 + dc*h;
    
    
    
    % Double Step
    c1a = c0 + dc*(h/2);
    kc = zeros(size(kinetic_comps));
    for i = 1:length(kc)
        kc(i) = c1a(ismember(components,kinetic_comps(i)));
    end
    vkg = evaluate_fluxes(kc,kinetic_growth);
    vkp = evaluate_fluxes(kc,kinetic_prod);
    vg = simulate_model(model,components,kinetic_comps,vkg,p_growth,f_growth,fix_flx,fix_comp);
    vp = simulate_model(model,components,kinetic_comps,vkp,p_prod,f_prod,fix_flx,fix_comp);
    
    f = phase_progress(c0,components,classification_data);
    v = (1-f)*vg + f*vp;
    dc = eval_reactor(c1a,v,c_in);
    c2 = c1a + dc*(h/2);
    

    % Check error
    rel_err = rms((c2-c1)/rtol);
    
    % Step Acceptance/Rejection
    if rel_err < 1
        % Accept Step
        ctr = ctr+1;
        c0 = c2;
        profiles(ctr,:) = c0;
        phase_transition(ctr) = f;
        flux_growth(ctr,:) = vg';
        flux_prod(ctr,:) = vp';
        t_curr = t_curr+h;
        time(ctr) = t_curr;
        h = h/max(sqrt(rel_err),0.5);
        h = min(h,time_range(2)-t_curr);
        msg = ['time = ',num2str(t_curr),' days of ',num2str(time_range(2)),' days.     f = ',num2str(f)];
        waitbar(t_curr/time_range(2),wbar,msg);
        
    else
        h = h/1.5;
    end
    
    if abs(time_range(2)-t_curr) < 1e-7
        done = true;
    end
end
delete(wbar);
dFBA_results.time = time(1:ctr);
dFBA_results.profiles = profiles(1:ctr,:);
dFBA_results.flux_growth = flux_growth(1:ctr,:);
dFBA_results.flux_prod = flux_prod(1:ctr,:);
dFBA_results.phase_transition = phase_transition(1:ctr);
dFBA_results.notes = {notes};
dFBA_results.condition = dFBA_data.Condition;

end

function v = simulate_model(model,rxns,rxn_in,v_in,priority,objs,fix_flx,fix_comp)
% Helper function to simulate the metabolic model using specified uptake
% rates and report secretion fluxes

% Update Reaction bounds
m1 = model;
for i = 1:length(rxn_in)
    m1.lb(ismember(m1.rxns,rxn_in(i))) = v_in(i);
    if ismember(rxn_in(i),fix_comp)
        m1.ub(ismember(m1.rxns,rxn_in(i))) = fix_flx(ismember(fix_comp,rxn_in(i)))*v_in(i);
    end
    %m1.ub(ismember(m1.rxns,rxn_in(i))) = 0.99*v_in(i);
end

for i = 1:length(priority)
    m1 = changeObjective(m1,priority(i),1);
    sol = optimizeCbModel(m1);
    m1.lb(ismember(m1.rxns,priority(i))) = objs(i)*sol.obj;
    %m1.ub(ismember(m1.rxns,priority(i))) = 1.01*m1.lb(ismember(m1.rxns,priority(i)));
end

sol = optimizeCbModel(m1,'min');
x = sol.x;

v = zeros(size(rxns));
for i = 1:length(v)
    v(i) = x(ismember(m1.rxns,rxns(i)));
end
v(3:end) = v(3:end)*0.35;

end

function f = phase_progress(c,comp_names,phase_classifier)
% Helper function to determine progress towards production phase.
% f is a quantity between 0 and 1 with 0 indicating all cells in growth
% phase and 1 indicating all cells in the production phase.
cx = zeros(size(phase_classifier.components));
for i = 1:length(cx)
    cx(i) = c(ismember(comp_names,phase_classifier.components(i)));
end
L = phase_classifier.Beta;
K = phase_classifier.Bias;
scale = phase_classifier.scaling;
L = L.*scale;
E = phase_classifier.log_coeff;

w = cx'*L + K;
w = w*E(2) + E(1);
f = exp(w)/(1+exp(w));
%f = 1/(1+exp(w));

end


function v = evaluate_fluxes(c,kinetic)
% Helper function to evaluate fluxes using concentrations and
% Michaelis-Menten Rate Laws
v = zeros(size(c));
for i = 1:length(v)
    v(i) = -kinetic(1,i)*c(i)/(kinetic(2,i)+c(i));
    v(i) = min(v(i),0);
end

end



function dc = eval_reactor(c,v,c_in)
% Helper function to evaluate ODE function (conservation of mass) 
D = 1;
c = max(c,0);
et = ones(size(c));
et(1) = 0;
dc = zeros(size(c));

for i = 1:numel(dc)
    if i < 3
        dc(i) = v(i)*c(i);
    else
        dc(i) = v(i)*c(1)+ D*(c_in(i)-(et(i)*c(i)));
    end
    
end
end









