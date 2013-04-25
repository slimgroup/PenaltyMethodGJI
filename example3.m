% This program is part of the paper
% "Mitigating local minima in full-waveform inversion by expanding the search space",
% T. van Leeuwen and F.J. Herrmann, 2013 (submitted to GJI).
%
% Copyright (C) 2013 Tristan van Leeuwen (tleeuwen@eos.ubc.ca)
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% define model
n = [101 301];
h = [50 50];
z = [0:n(1)-1]'*h(1);
x = [0:n(2)-1]*h(2);
e = ones(n(2),1);

v  = 2000 + 0.7*z + 200*exp(-1e-6*(z - 2000).^2);
m  = 1./v.^2;
v0 = 2000 + 0.7*z;
m0 = 1./v0.^2;

% inversion parameters
niter  = 1000;
L1     = 1e17;
L2     = 1e13;
lambda = 1e4;

% make data
model.n = n;
model.h = h;
model.f = 5;
model.Is = {2,1};
model.Ir = {2,1:301};
Q = 1;
D = F(m*e',Q,model);

% inversion with reduced misfit
m1 = m0;
for k = 1:niter
    [f1,g1] = misfit_red(m1*e',Q,D,model);
    info1(k,:) = [f1 norm(g1)];
    m1 = m1 - (1/L1)*reshape(g1,n)*e;
end

% inversion with penalty method
m2 = m0;
for k = 1:niter
    [f2,g2,h2] = misfit_pen(m2*e',Q,D,lambda,model);
    info2(k,:) = [f2 norm(g2)];
    m2 = m2 - (1/L2)*(reshape(g2,n)*e);
end

% data for final models
D1 = F(m1*e',Q,model);
D2 = F(m2*e',Q,model);

% plot results
figure;
plot(1./sqrt(m1)-v0,z,'b',1./sqrt(m2)-v0,z,'r',1./sqrt(m)-v0,z,'k--','linewidth',2);axis ij;xlim([-50 250]);set(gca,'plotboxaspectratio',[1 2 1],'fontsize',20);
xlabel('\delta v [m/s]','fontsize',20);ylabel('z [m]','fontsize',20);
legend('reduced','penalty');

figure;
loglog(1:niter,info1(:,1)/info1(1,1),'b',1:niter,info2(:,1)/info2(1,1),'r',1:niter,info1(:,2)/info1(1,2),'b--',1:niter,info2(:,2)/info2(1,2),'r--','linewidth',2);
axis tight;set(gca,'fontsize',20);
xlabel('iterations','fontsize',20);ylabel('magnitude','fontsize',20);
legend('reduced','penalty');

figure;
plot(model.Ir{2},real(D1),'b',model.Ir{2},real(D2),'r',model.Ir{2},real(D),'k--','linewidth',2);
set(gca,'fontsize',20,'plotboxaspectratio',[2 1 1]);
ylim([-1 1]);xlim([model.Ir{2}([1 end])]);
xlabel('receiver index','fontsize',20);ylabel('Re(d)','fontsize',20);
legend('reduced','penalty');

% save figures
print(1,'-depsc',['example3a']);
print(2,'-depsc',['example3b']);
print(3,'-depsc',['example3c']);
