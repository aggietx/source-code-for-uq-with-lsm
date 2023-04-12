% Test of feedback particle filter
rng(1); % for the sake of reproducing the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Generate the true signal %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nt = 5000;
dt = 0.005;
X_truth = zeros(1,Nt);
Z = zeros(1,Nt);
% Regime I: bimodal
a_X = 2.8; b_X = 8; c_X = 4; f_X = -0.4; sigma_X = 0.7;
% Regime II: nearly Gaussian
% a_X = 2.8; b_X = 4; c_X = 4; f_X = 0.4; sigma_X = 0.7;
sigma_Z = .1;
for i = 2:Nt
    hX = - Z(i-1) + X_truth(i-1)/5;
    X_truth(i) = X_truth(i-1) + (-a_X * X_truth(i-1) + b_X * X_truth(i-1)^2 - c_X * X_truth(i-1)^3 + f_X) * dt + sigma_X * sqrt(dt) * randn;
    Z(i) = Z(i-1) + hX * dt + sigma_Z * sqrt(dt) * randn;
end
figure
subplot(2,5,1:3)
plot(dt:dt:Nt*dt, X_truth,'b','linewidth',2)
box on
set(gca,'fontsize',12)
title('(a) Time series of state variable X')
subplot(2,5,4)
[fi,xx] = ksdensity(X_truth);
fi_G = normpdf(xx, mean(X_truth), std(X_truth));
hold on
plot(xx,fi,'b','linewidth',2)
plot(xx,fi_G,'--k','linewidth',2)
box on
set(gca,'fontsize',12)
title('(b) PDF')
subplot(2,5,5)
hold on
plot(xx,fi,'b','linewidth',2)
plot(xx,fi_G,'--k','linewidth',2)
set(gca,'yscale','log')
box on
set(gca,'fontsize',12)
legend('Truth','Gaussian fit')
ylim([1e-3,5])
title('(c) PDF in log scale')

subplot(2,5,6:8)
plot(dt:dt:Nt*dt, Z,'b','linewidth',2)
box on
set(gca,'fontsize',12)
title('(d) Time series of obs Z')
subplot(2,5,9)
[fi,xx] = ksdensity(Z);
fi_G = normpdf(xx, mean(Z), std(Z));
hold on
plot(xx,fi,'b','linewidth',2)
plot(xx,fi_G,'--k','linewidth',2)
box on
set(gca,'fontsize',12)
title('(e) PDF')
subplot(2,5,10)
hold on
plot(xx,fi,'b','linewidth',2)
plot(xx,fi_G,'--k','linewidth',2)
set(gca,'yscale','log')
box on
set(gca,'fontsize',12)
% legend('Truth','Gaussian fit')
ylim([1e-3,5])
title('(f) PDF in log scale')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Feedback particle filter %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 100; L = 30; epsilon = 0.05;
X_FPF_all = zeros(N,Nt);
K_FPF_all = zeros(N,Nt);
X = X_truth(1) + randn(N,1) * 0.1; % initial ensembles
X_FPF_all(:,1) = X;
K = zeros(N,1);
Phi_prev = ones(N,1);
for t = 2: Nt
    hX = - Z(t-1) + X/5; 
    g = exp(- (ones(N,1)*X' - X*ones(1,N)).^2/4/epsilon);
    k = g./(sqrt(sum(g,2) * sum(g,1)));
    T = k./ (sum(k,2) * ones(1,N));
    h_hat = mean(hX);
    Phi = Phi_prev;
    for l = 1:L
        Phi = T * Phi + epsilon * (hX - h_hat);
        Phi = Phi - mean(Phi);
    end
    r = Phi + epsilon * (hX - h_hat);
    a = 1/2/epsilon * T .* (ones(N,1)*r' - T*r*ones(1,N));
    K = a * X;        
    if kurtosis(K)>30
        K = mean(K) * ones(N,1);
    end

    X_new = X + (-a_X * X + b_X * X.^2 - c_X * X.^3 + f_X) * dt + sigma_X * sqrt(dt) * randn(N,1) ...
        + K / sigma_Z^2 .* ((Z(t) - Z(t-1))/dt - 1/2 * (hX + h_hat)) * dt;
    X = X_new;
    if mod(t,10 ) == 1

    end
    Phi_prev = Phi;
    if t < Nt
        X_FPF_all(:,t) = X;
    end
    K_FPF_all(:,t) = K;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Ensemble Kalman-Bucy filter %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X_EnKBF_all = zeros(N,Nt);
K_EnKBF_all = zeros(1,Nt);
X = X_truth(1) + randn(N,1) * 0.1; % initial ensembles
X_EnKBF_all(:,1) = X;
for t = 2: Nt
    hX = -Z(t-1) + X/5; 
    K = sum((X - mean(X)) .* (hX - mean(hX)))/(N-1);
    X = X + (-a_X * X + b_X * X.^2 - c_X * X.^3 + f_X) * dt + sigma_X * sqrt(dt) * randn(N,1) ...
        +  K / sigma_Z^2 .* (Z(t) - Z(t-1) - hX * dt - randn(N,1) * sigma_Z * sqrt(dt));
    if t < Nt
        X_EnKBF_all(:,t) = X;
    end
    K_EnKBF_all(t) = K;
end
% posterior mean trajectory
figure
hold on
plot(dt:dt:Nt*dt, X_truth, 'b', 'linewidth', 2);
plot(dt:dt:Nt*dt, mean(X_FPF_all), 'r', 'linewidth', 2);
plot(dt:dt:Nt*dt, mean(X_EnKBF_all), 'g', 'linewidth', 2);
box on
set(gca,'fontsize',12)
legend('Truth','FPF','EnKBF')
title('Truth and posterior mean time series')
% posterior distribution
figure
for i = 1:4
    subplot(2,2,i)
    hold on
    [fi_FPF,xx] = ksdensity(X_FPF_all(:,100*i));    
    [fi_EnKBF,xx] = ksdensity(X_EnKBF_all(:,100*i),xx);
    temp = max(fi_FPF);
    plot([X_truth(100*i),X_truth(100*i)],[0,temp],'b','linewidth',2)
    plot(xx,fi_FPF,'r','linewidth',2)
    plot(xx,fi_EnKBF,'g','linewidth',2)    
    if i == 1
        legend('Truth','FPF','EnKBF')
    end
    box on
    set(gca,'fontsize',12)
    title(['PDF at t = ', num2str(round(i*100*dt))])
end
% 'Kalman gain'
figure
for i = 1:4
    subplot(2,2,i)
    hold on
    [fi_FPF,xx] = ksdensity(K_FPF_all(:,100*i));   
    temp = max(fi_FPF);
    plot([K_EnKBF_all(100*i),K_EnKBF_all(100*i)],[0,temp],'g','linewidth',2)
    plot(xx,fi_FPF,'r','linewidth',2)  
    if i == 1
        legend('EnKBF','FPF')
    end
    box on
    set(gca,'fontsize',12)
    title(['PDF at t = ', num2str(round(i*100*dt))])
end