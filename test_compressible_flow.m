rng(1);
n = 8;24; % model dim
L = 10; % num of obs
epsilon = 1;
k = [1, 0, 1, -1;%, 2, 0, 1, 2, 2, -1, 2, -2;
     0, 1, 1,  1];%, 0, 2, 2, 1, -1, 2, 2,  2];
kk = [k,k,k];

sigmal = 0.2;
omegak_p = 1/epsilon * sqrt(k(1,:).^2 + k(2,:).^2 + 1);
omegak_m = - 1/epsilon * sqrt(k(1,:).^2 + k(2,:).^2 + 1);

rk1 = [1./sqrt(k(1,:).^2 + k(2,:).^2+1) .* (-1i * k(2,:));
    1./sqrt(k(1,:).^2 + k(2,:).^2+1) .* (1i * k(1,:))];
rk2 = [1./sqrt(k(1,:).^2 + k(2,:).^2)/sqrt(2)/sqrt(k(1,:).^2 + k(2,:).^2 + 1) .* (1i * k(2,:) + k(1,:) .* sqrt(k(1,:).^2 + k(2,:).^2 + 1));
    1./sqrt(k(1,:).^2 + k(2,:).^2)/sqrt(2)/sqrt(k(1,:).^2 + k(2,:).^2 + 1) .* (-1i * k(1,:) + k(2,:) .* sqrt(k(1,:).^2 + k(2,:).^2 + 1))];
rk3 = [1./sqrt(k(1,:).^2 + k(2,:).^2)/sqrt(2)/sqrt(k(1,:).^2 + k(2,:).^2 + 1) .* (1i * k(2,:) - k(1,:) .* sqrt(k(1,:).^2 + k(2,:).^2 + 1));
    1./sqrt(k(1,:).^2 + k(2,:).^2)/sqrt(2)/sqrt(k(1,:).^2 + k(2,:).^2 + 1) .* (-1i * k(1,:) - k(2,:) .* sqrt(k(1,:).^2 + k(2,:).^2 + 1))];
rk = [rk1,rk2,rk3];
dk1 = -[0.5:0.05:1.05];
dk2 = -[0.65:0.05:1.2];
dk3 = -[0.65:0.05:1.2];
ak1 = zeros(n,n);
ak2 = zeros(n,n);
ak3 = zeros(n,n);
for i = 1:2:n
    ak1(i,i) = dk1((i+1)/2);
    ak1(i+1,i+1) = ak1(i,i);    
    ak2(i,i) = dk2((i+1)/2);
    ak2(i+1,i+1) = ak2(i,i); 
    ak2(i,i+1) = -omegak_p((i+1)/2);
    ak2(i+1,i) =  omegak_p((i+1)/2);
    ak3(i,i) = dk3((i+1)/2);
    ak3(i+1,i+1) = ak3(i,i);
    ak3(i,i+1) = -omegak_m((i+1)/2);
    ak3(i+1,i) =  omegak_m((i+1)/2);
end


a0 = zeros(3*n,1);
ak = blkdiag(ak1,ak2,ak3);
dt = 0.0001;
a1 = ak;
sigmak1 = eye(n);
sigmak2 = eye(n);
sigmak3 = eye(n);
sigmak = blkdiag(sigmak1,sigmak2,sigmak3);

b2 = sigmak;
Sgm1 = sigmal * eye(2*L);
InvBoB = (Sgm1 * Sgm1')^(-1);

T = 5;
NN = round(T/dt);
uII = zeros(3*n,NN);
uI = zeros(2*L,NN);
G = zeros(L*2,3*n);
uI(:,1) = rand(2*L,1)*pi*2-pi;
x = uI(1:2:end,1);x=x';
y = uI(2:2:end,1);y=y';
rng(2)
for i = 2:NN
    uII(:,i) = uII(:,i-1) + (a0 + a1 * uII(:,i-1)) * dt + b2 * sqrt(dt) * randn(3*n,1); % true signal up to current
end


for i = 2:NN
    G(1:2:end,1:2:end) = real(exp(1i*(x'*kk(1,:)+y'*kk(2,:))).*(ones(L,1) * rk(1,:)));
    G(1:2:end,2:2:end) = real(exp(1i*(x'*kk(1,:)+y'*kk(2,:))).* 1i.*(ones(L,1) *rk(1,:)));
    G(2:2:end,1:2:end) = real(exp(1i*(x'*kk(1,:)+y'*kk(2,:))).*(ones(L,1) *rk(2,:)));
    G(2:2:end,2:2:end) = real(exp(1i*(x'*kk(1,:)+y'*kk(2,:))).* 1i.*(ones(L,1) *rk(2,:)));
    G = G*2;
    A1 = G;
    uI(:,i) = uI(:,i-1) + A1 * uII(:,i-1) * dt + Sgm1 * randn(2*L,1) * sqrt(dt); % obs up to current
    x = uI(1:2:end,i-1);x=x';
    y = uI(2:2:end,i-1);y=y';
end

Rm = -diag(diag(sigmak).^2/2./diag(ak));

mu0 = zeros(3*n,1);
R0 = eye(3*n);
u_post = zeros(3*n,NN);
u_post(:,1) = mu0;
u_cov = zeros(3*n,NN);
u_cov(:,1) = diag(R0);
info_gain_total = zeros(1,NN);
info_gain_signal = zeros(1,NN);
info_gain_dispersion = zeros(1,NN);
for i = 2:NN
    x = uI(1:2:end,i-1);x=x';
    y = uI(2:2:end,i-1);y=y';
    G(1:2:end,1:2:end) = real(exp(1i*(x'*kk(1,:)+y'*kk(2,:))).*(ones(L,1) * rk(1,:)));
    G(1:2:end,2:2:end) = real(exp(1i*(x'*kk(1,:)+y'*kk(2,:))).* 1i.*(ones(L,1) *rk(1,:)));
    G(2:2:end,1:2:end) = real(exp(1i*(x'*kk(1,:)+y'*kk(2,:))).*(ones(L,1) *rk(2,:)));
    G(2:2:end,2:2:end) = real(exp(1i*(x'*kk(1,:)+y'*kk(2,:))).* 1i.*(ones(L,1) *rk(2,:)));
    G = G*2;
    A1 = G;
    mu = mu0 + (a0 + a1 * mu0) * dt + (R0 * A1') * InvBoB * (uI(:,i)-uI(:,i-1) - A1 * mu0 * dt);
    R = R0 + (a1 * R0 + R0* a1' + b2*b2' - (R0*A1') * InvBoB * (R0*A1')')*dt;
    u_post(:,i) = mu;
    u_cov(:,i) = diag(R);
    mu0 = mu;
    R0 = R;

    info_gain_signal(i) = 1/2*(mu)'*Rm^(-1)*(mu);
    info_gain_dispersion(i) =  -1/2 * log(det(R/Rm)) + 1/2 * trace(R/Rm) - 1/2 * 3*n;
    info_gain_total(i) = info_gain_signal(i) + info_gain_dispersion(i);
end
yval_current = info_gain_total(end);
numtemp = [1, 3, 5, 11];
figure
for i = 1:3
    for k = 1:4
        subplot(4,3,(k-1)*3+i)
        hold on
        plot(dt:dt:NN*dt,uII(numtemp(k)*2-1 + (i-1)*24,:),'b')
        plot(dt:dt:NN*dt,u_post(numtemp(k)*2-1 + (i-1)*24,:),'r')
        box on
        set(gca,'fontsize',12)
        ylim([-2,2])
        if i == 1
            if k == 1
                text(-2,0,'mode (1,0)','fontsize',14)
                legend('Truth','Assimilated mean')
            elseif k == 2
                text(-2,0,'mode (1,1)','fontsize',14)
            elseif k == 3
                text(-2,0,'mode (2,0)','fontsize',14)
            elseif k == 4
                text(-2,0,'mode (2,2)','fontsize',14)
            end
        end
        if k == 1
            if i == 1
                title('(a) GB modes')
            elseif i == 2
                title('(b) Gravity modes +')
            else
                title('(c) Gravity modes -')
            end
        end
        if k == 4
            xlabel('t')
        end
    end
end


uI_new0 = uI(:,end);
uII_new0 = uII(:,end);
mu_new0 = mu;
R_new0 = R;



% 

nb = 20;
uI = mod(uI+pi,2*pi)-pi;
[xx,yy] = meshgrid(linspace(-pi,pi,nb),linspace(-pi,pi,nb));
v1 = zeros(nb,nb);
v2 = zeros(nb,nb);
v1f = zeros(nb,nb);
v2f = zeros(nb,nb);
GG = zeros(2,n*3);
figure
for i = 1:3
    subplot(2,3,i)
    time = i * 20000-10000;
    u_all = [uII(1:2:end,time);uII(2:2:end,time)];
    uf_all = [u_post(1:2:end,time);u_post(2:2:end,time)];
    for j1 = 1:nb
        for j2 = 1:nb
            x = xx(j1,j2);
            y = yy(j1,j2);
            GG(1:2:end,1:2:end) = real(exp(1i*(x'*kk(1,:)+y'*kk(2,:))).*(rk(1,:)));
            GG(1:2:end,2:2:end) = real(exp(1i*(x'*kk(1,:)+y'*kk(2,:))).* 1i.*(rk(1,:)));
            GG(2:2:end,1:2:end) = real(exp(1i*(x'*kk(1,:)+y'*kk(2,:))).*(rk(2,:)));
            GG(2:2:end,2:2:end) = real(exp(1i*(x'*kk(1,:)+y'*kk(2,:))).* 1i.*(rk(2,:)));
            GG = GG*2;
            A1G = GG;
            v = A1G * u_all;
            vf = A1G * uf_all;
            v1(j1,j2) = v(1);
            v2(j1,j2) = v(2);
            v1f(j1,j2) = vf(1);
            v2f(j1,j2) = vf(2);
        end
    end
    
    v_amp = sqrt(v1.^2 + v2.^2);
    hold on
    contourf(xx,yy,v_amp,'linestyle','none')
    quiver(xx,yy,v1,v2,'k','linewidth',2)
    plot(uI(1:2:end,time),uI(2:2:end,time),'mo','linewidth',3)
    xlim([-pi,pi]);
    ylim([-pi,pi]);
    title(['Time = ',num2str(time * dt)],'fontsize',14);
    set(gca,'fontsize',12)
    box on
    xlabel('x','fontsize',12)
    ylabel('y','fontsize',12)
%     colormap parula
    
    subplot(2,3,i+3)
%     time = i * 10000;
    
    vf_amp = sqrt(v1f.^2 + v2f.^2);
    hold on
    contourf(xx,yy,vf_amp,'linestyle','none')
    quiver(xx,yy,v1f,v2f,'k','linewidth',2)
    plot(uI(1:2:end,time),uI(2:2:end,time),'mo','linewidth',3)
    xlim([-pi,pi]);
    ylim([-pi,pi]);
    title(['Time = ',num2str(time * dt)],'fontsize',14);
    set(gca,'fontsize',12)
    box on
    xlabel('x','fontsize',12)
    ylabel('y','fontsize',12)
%     colormap winter
end




