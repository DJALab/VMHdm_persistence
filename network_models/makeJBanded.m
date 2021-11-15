function J = makeJBanded(N,sp,sig,offset)

J = zeros(N);
for n=1:N
    band = normpdf((1:N),n+offset,sig);
    band = band/max(band);
    J(n,:)   = rand(1,N).*(rand(1,N)<(sp/mean(band).*band))/sqrt(N*sp);
end