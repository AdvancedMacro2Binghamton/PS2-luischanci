%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Binghamton University                         %
% PhD in Economics                              %
% ECON634 Advanced Macroeconomics               %
% P.S. 2                                        %
% Fall 2017                                     %
% Luis Chancí (lchanci1@binghamton.edu)         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; a  = 0.35; b = 0.99; d = 0.025; s = 2;
z  = [.678 1.1]';                        % A^l and A^h
p  = [.926 (1-.977)]';                   % pi(ll) and pi(hl)
kss= (a*inv(inv(b)-(1-d)))^inv(1-a);     % I use the k(steady state)
k  = linspace(0.001*kss, 1.3*kss, 1000); %   between 0.1% and 130% of kss
km = repmat(k',[1 length(k)]);
c  = cat(3,z(1)*km.^a+(1-d)*km-km',z(2)*km.^a+(1-d)*km-km'); 
u  = c.^(1-s)/(1-s); u(c<0)=-Inf;
v0 = zeros(length(z),length(k));
e  = 1; t0 = tic;
while e>1e-06
    for m=1:length(z)
        v(:,:,m) = u(:,:,m)+b*(p(m)*repmat(v0(1,:), [length(k) 1])+...
                           (1-p(m))*repmat(v0(2,:), [length(k) 1]));
        [vfn(m,:),idx(m,:)] = max(v(:,:,m),[],2);
    end
    e = max(max(abs(vfn - v0)));
    v0 = vfn;
end; t = toc(t0); g = k(idx);   % Now we can plot

% Simulation.
S = bsxfun(@gt,[.926;.977],rand(1,length(k)));
y = [];  j = 1;
for i = 1:length(k)
    y = [y z(j)*g(j,i)^a];
    if S(j,i) == 1
         if j == 1; j = 2; elseif j == 2; j = 1; end
    end
end                       % Now we can get var(y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%