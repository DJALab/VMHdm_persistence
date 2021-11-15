function J = decorrelateWeights(J,g,dt,tauS,tau,thr,beta,tmax,scaleEvery)

N = size(J,1);
doPlot=0;
if(doPlot)
    figure(123);clf;
    subplot(2,1,1);
    h=image(J(1:100,1:100)*640);
    h2=title('0');
    colorbar;
    subplot(2,1,2);
    h3=plot(1:N,1:N);
    ylim([0 20]);
end
rsm = zeros(N,1);
[r,p,x] = deal(zeros(N,1));
for t=1:tmax
    p   = r + (p - r) * exp(-dt/tauS);
    Iin = max(g*J*p + 0.085 + rand(N,1)*5*thr/sqrt(t),0); %drive the network with random input
    x   = Iin + (x - Iin) * exp(-dt/tau);

    r  = (x >= thr)*1/dt/100;
    x(r~=0) = 0;
    dJ = zeros(N);
    for n=find(r~=0)'
        dJ(n,:) = -rsm';
        dJ(n,rsm==0) = sum(rsm)/sum(rsm==0);
    end
    rsm = rsm*exp(-dt/beta) + r;
    J = J + dJ*.0001;
    J = max(J-diag(diag(J)),0);
    if(doPlot)
        set(h,'cdata',J(1:100,1:100)*640);
        set(h2,'string',t);
        set(h3,'ydata',rsm);
        drawnow;
    end
    if(mod(t,scaleEvery)==0)
        [e,v]=eig(J-mean(J(:)));
        J=J/max(real(diag(v)));
    end
end

[e,v]=eig(J-mean(J(:)));
J=J/max(real(diag(v)));