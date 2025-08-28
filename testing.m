%% Instantiate parameter structure

parameters = struct('s1',0, ...     % beginning of signalling region
                    's2',0.3, ...   % end of signalling region
                    'r1',0.4, ...   % beginning of responsive region
                    'r2',0.5, ...   % end of responsive region
                    'dt',0.001, ... % simulation time step
                    'tfin',100, ... % ending time
                    'del',0.1, ...  % sharpness of hyperbolic cutoffs
                    'eps',0.025, ...% diffusive coefficient
                    'alph',0.9, ... % linear inhibition term
                    'fr',2, ...     % frame rate
                    'pr',1);        % pause rate

%%

parameters.s1    = 0.1;
parameters.s2    = 0.3;
parameters.r1    = 0.4;
parameters.r2    = 0.5;
parameters.dt    = 0.001;
parameters.tfin  = 0.5;
parameters.del   = 0.1;
parameters.eps   = 0.025;
parameters.alph  = 0.9;
parameters.fr    = 100; 
parameters.pr    = 0.01;

disp('start')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test 1

n = 10;
x = linspace(0,1,2^n);

% test fatten_points_polynomial
d_0 = fatten_points_polynomial(x,0.5,0.3);
d_0 = d_0.';

% test pseudospectral
[ps_t, ps_d] = pseudospectral(d_0,parameters);

% test process_pseudospectral
ps_f = process_pseudospectral(x,ps_t,ps_d);

runs = 2;
pts = 200;

rr = zeros([2*runs,pts]);

e_0 = zeros([2*runs,pts]);

for k = 1:runs

    % test construct_em
    p = construct_em(x,d_0,pts);
    % test forward_euler_noise
    [en_t, en_d] = forward_euler_noise(p,parameters);
    % test forward_euler
    [e_t, e_d] = forward_euler(p,parameters);

    e_0(k,:) = en_d(1,:);
    e_0(2+k,:) = e_d(1,:);
    rr(k,:) = en_d(1,end,:);

end

mms = zeros([runs,3]);
mms_0 = zeros([runs,3]);
ff_d = zeros([runs,1,length(x)]);
ic_f = zeros([1+runs,1,length(x)]);
ic_f(1,1,:) = d_0;
for k = 1:runs

    % test metric_lp_1
    mms(k,:) = metric_lp_1(rr(k,:),x,ps_d(end,:),1);
    ho = fatten_points_polynomial(x,rr(k,:),mms(k,2));
    ff_d(k,1,:) = ho;
    
    mms_0(k,:) = metric_lp_1(e_0(k,:),x,d_0.',1);
    ho = fatten_points_polynomial(x,e_0(k,:),mms_0(k,2));
    ic_f(k+1,1,:) = ho;

end

al_f = cat(1,ps_f(:,end,:),ff_d);
names = ["ps","en1","en2","e3","e4"];
colors = [linspace(0,1,runs);zeros([1,runs]);linspace(1,0,runs)].';
frame(x,al_f,parameters,2,'Names',names,'Colors',colors,'Time',parameters.tfin);
frame(x,ic_f,parameters,3,'Names',names,'Colors',colors,'Time',0);

disp('Done Test 1')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test 2

particles = rand([1, 100]);
n = 10;
x = linspace(0,1,2^n);

particles_fat = fatten_points_polynomial(x,particles,0.1);

[ps_t, ps_d] = pseudospectral(particles_fat.',parameters);
ps_f = process_pseudospectral(x,ps_t,ps_d);

[en_t, en_d] = forward_euler_noise(particles,parameters);
en_f = process_euler(x,en_t,en_d,0.1);

[fe_t, fe_d] = forward_euler(particles,parameters);
fe_f = process_euler(x,fe_t,fe_d,0.1);

al_f = cat(1,fe_f,ps_f);
al_f = cat(1,al_f,en_f);
names = ["forward euler","pseudospectral","noisy euler"];
animate(x,al_f,parameters,4,'Names',names);

disp('Done Test 2')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test 3

n = 10;
x = linspace(0,1,2^n);

d_0 = fatten_points_polynomial(x,0.5,0.3)+ 0.2;
d_0 = d_0/lp_integrate(x,d_0,1);
d_0 = d_0.';

e_0 = construct_em(x,d_0,1000);
e_0_t = fatten_points_polynomial(x,e_0,0.001)/10;
e_t = fatten_points_polynomial(x,e_0,0.1);

figure(5)
clf
hold on

plot(x,d_0);
plot(x,e_0_t)
plot(x,e_t)

disp('Done Test 3')

disp('end')