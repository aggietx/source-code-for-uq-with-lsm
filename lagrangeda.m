clear
close all;
xMin = 0;
xMax = 2*pi;

yMin = 0;
yMax = 2*pi;

nGridpoints = 50;

x = linspace(0, 2*pi, nGridpoints);
y = linspace(0, 2*pi, nGridpoints);
[X, Y] = meshgrid(x, y);

maxWavenumber = 15;
wavenumbers = -maxWavenumber:maxWavenumber;
nWavenumbers = length(wavenumbers);

basisFunctions = nan(nWavenumbers, nWavenumbers, nGridpoints, nGridpoints, 2);

for iWavenumber = 1:nWavenumbers
    for jWavenumber = 1:nWavenumbers
        % Notation comes from Mohamad, Majda (2020)
        
        k = [wavenumbers(iWavenumber); wavenumbers(jWavenumber)];
        
        if norm(k) == 0
            r_k = [0; 0];
        else
            k_perp = [-k(2); k(1)];
            r_k = 1i*k_perp/norm(k);
        end
        
        Psi_k1 = r_k(1)*exp(1i*(k(1)*X + k(2)*Y));
        Psi_k2 = r_k(2)*exp(1i*(k(1)*X + k(2)*Y));
        
        basisFunctions(iWavenumber, jWavenumber, :, :, 1) = Psi_k1;
        basisFunctions(iWavenumber, jWavenumber, :, :, 2) = Psi_k2;
    end
end

alpha = 3;
E_0 = 1;
k_0 = 2;
energies = nan(nWavenumbers, nWavenumbers);
for iWavenumber = 1:nWavenumbers
    for jWavenumber = 1:nWavenumbers    
        k = [wavenumbers(iWavenumber); wavenumbers(jWavenumber)];
        if norm(k) <= k_0
            energies(iWavenumber, jWavenumber) = norm(k)*E_0;
        else
            energies(iWavenumber, jWavenumber) = k_0*E_0*abs(norm(k)/k_0)^-alpha;
        end
    end    
end










% OU Process
dampingCoeffs = zeros(nWavenumbers, nWavenumbers);
forcings = zeros(nWavenumbers, nWavenumbers);
noiseCoeffs = zeros(nWavenumbers, nWavenumbers);
means = zeros(nWavenumbers, nWavenumbers);
variances = zeros(nWavenumbers, nWavenumbers);

for iWavenumber = 1:nWavenumbers
    for jWavenumber = 1:nWavenumbers
        d = 1;
        nu = 1;
        k = [wavenumbers(iWavenumber); wavenumbers(jWavenumber)];
        d_k = d + nu*norm(k)^2;
        dampingCoeffs(iWavenumber, jWavenumber) = d_k;

        f_k = 0;
        forcings(iWavenumber, jWavenumber) = f_k;

        E_k = energies(iWavenumber, jWavenumber);
        sigma_k = sqrt(2*(2*E_k-f_k/d_k^2)*d_k);
        %sigma_k = sqrt(2*E_k*d_k);

        noiseCoeffs(iWavenumber, jWavenumber) = sigma_k;

        means(iWavenumber, jWavenumber) = f_k/d_k;
        variances(iWavenumber, jWavenumber) = sigma_k^2/(2*d_k);
    end
end

amplitudes = nan(nWavenumbers, nWavenumbers);
for iWavenumber = 1:nWavenumbers
    for jWavenumber = 1:nWavenumbers
        vbar = means(iWavenumber, jWavenumber);
        vsigma = variances(iWavenumber, jWavenumber);
        v_k = vbar + vsigma*(randn(1) + 1i*randn(1))/sqrt(2);
        amplitudes(iWavenumber, jWavenumber) = v_k;
    end    
end

for iWavenumber = 1:nWavenumbers
    for jWavenumber = 1:nWavenumbers
        amplitudes(nWavenumbers-iWavenumber+1, nWavenumbers-jWavenumber+1) = conj(amplitudes(iWavenumber, jWavenumber));
    end    
end

nTracers = 20;
positions = nan(2, nTracers);
for iTracer = 1:nTracers
    positions(:, iTracer) = rand(2, 1)*(xMax-xMin) + xMin;
end


dt = 0.001;
maxTime = 1;
times = 0:dt:maxTime;
nTimes = length(times);
nSkip = 0.01/dt;

sigma_x = 0.01;

figure
for iTime = 1:nTimes
    for iWavenumber = 1:nWavenumbers
        for jWavenumber = 1:nWavenumbers
            d_k = dampingCoeffs(iWavenumber, jWavenumber);
            f_k = forcings(iWavenumber, jWavenumber);
            sigma_k = noiseCoeffs(iWavenumber, jWavenumber);
            v_k = amplitudes(iWavenumber, jWavenumber);
            amplitudes(iWavenumber, jWavenumber) = v_k + (-d_k*v_k + f_k)*dt + sigma_k*sqrt(dt)*(randn(1) + 1i*randn(1))/sqrt(2);
        end    
    end
    
    for iWavenumber = 1:nWavenumbers
        for jWavenumber = 1:nWavenumbers
            amplitudes(nWavenumbers-iWavenumber+1, nWavenumbers-jWavenumber+1) = conj(amplitudes(iWavenumber, jWavenumber));
        end    
    end
    
    U = zeros(nGridpoints, nGridpoints);
    V = zeros(nGridpoints, nGridpoints);
    for iWavenumber = 1:nWavenumbers
        for jWavenumber = 1:nWavenumbers
            v_k = amplitudes(iWavenumber, jWavenumber);
            Psi_k1 = squeeze(basisFunctions(iWavenumber, jWavenumber, :, :, 1));
            Psi_k2 = squeeze(basisFunctions(iWavenumber, jWavenumber, :, :, 2));
            U = U + v_k*Psi_k1;
            V = V + v_k*Psi_k2;
        end    
    end

    u = interp2(X, Y, real(U), positions(1, :), positions(2, :));
    v = interp2(X, Y, real(V), positions(1, :), positions(2, :));
    velocities = [u; v];
    
    for iTracer = 1:nTracers
        positions(:, iTracer) = positions(:, iTracer) + velocities(:, iTracer)*dt+sigma_x*randn(2, 1)*sqrt(dt);
    end

    positions(1, :) = mod(positions(1, :) - xMin, xMax) + xMin;
    positions(2, :) = mod(positions(2, :) - yMin, yMax) + yMin;

%     if mod(iTime, nSkip) == 0 
    if iTime== nSkip    
        clf
        hold on
        quiver(X, Y, real(U), real(V))
        axis([xMin, xMax, yMin, yMax])
        s = scatter(positions(1, :), positions(2, :), 'r.');
        s.SizeData = 100;
        title(['t = ', num2str(times(iTime))])
        pbaspect([1, 1, 1])
        pause(0.1)
    end
end

energy = [];
trueEnergy = [];
radii = [];
for iWavenumber = 1:nWavenumbers
    for jWavenumber = 1:nWavenumbers
        k = [wavenumbers(iWavenumber); wavenumbers(jWavenumber)];
        radii = [radii; norm(k)];
        energy = [energy; 0.5*abs(amplitudes(iWavenumber, jWavenumber))^2];
        trueEnergy = [trueEnergy, energies(iWavenumber, jWavenumber)];
    end
end

[radii, idx] = unique(radii);
energy = energy(idx);
trueEnergy = trueEnergy(idx);

figure
loglog(radii, energy, radii, trueEnergy);%, radii, A)
%A = energy;

% OU Process









