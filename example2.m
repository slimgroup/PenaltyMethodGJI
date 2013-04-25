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
n = [101 201];
h = [50 50];
[zz,xx] = ndgrid([0:n(1)-1]*h(1),[0:n(2)-1]*h(2));

v = @(v0,alpha)(v0 + alpha*zz(:));

% make data
model.n = n;
model.h = h;
model.f = 5;
model.Is = {2,1};
model.Ir = {2,1:201};
Q = 1;
D = F(1./v(2000,0.75).^2,Q,model);

% scan over v0
vs = [1750:25:2250];

for k = 1:length(vs);
    vk = v(vs(k),0.75);
    
    fv_red(k)  = misfit_red(1./vk.^2,Q,D,model);
    fv_pen1(k) = misfit_pen(1./vk.^2,Q,D,1e3,model);
    fv_pen2(k) = misfit_pen(1./vk.^2,Q,D,1e6,model);
    fv_pen3(k) = misfit_pen(1./vk.^2,Q,D,1e9,model);
end
   
% scan over alpha
as = [0.5:.02:1];

for k = 1:length(as);
    vk = v(2000,as(k));
    
    fa_red(k)  = misfit_red(1./vk.^2,Q,D,model);
    fa_pen1(k) = misfit_pen(1./vk.^2,Q,D,1e3,model);
    fa_pen2(k) = misfit_pen(1./vk.^2,Q,D,1e6,model);
    fa_pen3(k) = misfit_pen(1./vk.^2,Q,D,1e9,model);
end

% plot results
figure;plot(vs,fv_pen1/max(fv_pen1),'r',vs,fv_pen2/max(fv_pen2),'g',vs,fv_pen3/max(fv_pen3),'b',vs,fv_red/max(fv_red),'k--','linewidth',2);
xlabel('v0 [m/s]','fontsize',20);ylabel('misfit','fontsize',20);legend('\lambda = 10^3','\lambda = 10^6','\lambda = 10^9','reduced');
set(gca,'fontsize',20);axis tight

figure;plot(as,fa_pen1/max(fa_pen1),'r',as,fa_pen2/max(fa_pen2),'g',as,fa_pen3/max(fa_pen3),'b',as,fa_red/max(fa_red),'k--','linewidth',2);
xlabel('\alpha [1/s]','fontsize',20);ylabel('misfit','fontsize',20);legend('\lambda = 10^3','\lambda = 10^6','\lambda = 10^9','reduced');
set(gca,'fontsize',20);axis tight

% save figures
print(1,'-depsc',['example2a']);
print(2,'-depsc',['example2b']);