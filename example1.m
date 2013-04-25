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

% define velocity model
n  = [51 51];
h  = [20 20];
z  = [0:n(1)-1]'*h(1);
x  = [0:n(2)-1]*h(2);

v0 = 2000*ones(n);
dv = zeros(n);
dv(21:31,21:31) = 100;

m  = 1./(v0(:) + dv(:)).^2;
m0 = 1./(v0(:)).^2;

% frequency
f  =  10;

% sampling operators
Ps = getP(n,[1:51],2);
Pr = getP(n,[1:51],50);

% get Helmholtz matrix and source functions
A  = getA(f,m,h,n);
Q  = speye(51);

% compute wavefield and data
U  = A\Ps'*Q;
D  = Pr*U;

% set penalty parameter
lambda = 1;

% get Helmholtz operator for constant background model, m0
A0 = getA(f,m0,h,n);

% compute wavefield for m0
U0 = A0\Ps'*Q;

% compute wavefield from data-augmented wave-equation
U1 = [lambda*A0;Pr]\[lambda*Ps'*Q;D];

% solve for model from equation (7)
m1 = m0 - sum(conj(U1).*(A0*U1 - Ps'*Q),2)./((2*pi*f)^2*sum(abs(U1).^2,2));

% repeat procedure 4 more times, includes windowing operator
mk = m1;
w  = ones(n); w(1:51,50) = 0;
W  = spdiags(w(:),0,prod(n),prod(n));
for k = 1:9
    Ak = getA(f,mk,h,n);
    Uk = [lambda*Ak;Pr]\[lambda*Ps'*Q;D];
    mk = mk - W*(sum(conj(Uk).*(Ak*Uk - Ps'*Q),2)./((2*pi*f)^2*sum(abs(Uk).^2,2)));
end

% plotting
figure;imagesc(x,z,v0+dv,[2000 2100]);colorbar;
xlabel('x [m]','fontsize',20);ylabel('z [m]','fontsize',20);
set(gca,'fontsize',20);

figure;imagesc(x,z,reshape(real(U(:,26)-U0(:,26)),n));colormap(gray)
xlabel('x [m]','fontsize',20);ylabel('z [m]','fontsize',20);
set(gca,'fontsize',20);

figure;imagesc(x,z,reshape(real(U1(:,26)-U0(:,26)),n));colormap(gray)
xlabel('x [m]','fontsize',20);ylabel('z [m]','fontsize',20);
set(gca,'fontsize',20);

figure;imagesc(x,z,reshape(real(1./sqrt(m1)),n),[2000 2100]);colorbar;
xlabel('x [m]','fontsize',20);ylabel('z [m]','fontsize',20);
set(gca,'fontsize',20);

figure;imagesc(x,z,reshape(real(Uk(:,26)-U0(:,26)),n));colormap(gray)
xlabel('x [m]','fontsize',20);ylabel('z [m]','fontsize',20);
set(gca,'fontsize',20);

figure;imagesc(x,z,reshape(real(1./sqrt(mk)),n),[2000 2100]);colorbar;
xlabel('x [m]','fontsize',20);ylabel('z [m]','fontsize',20);
set(gca,'fontsize',20);


% save figures
print(1,'-depsc',['example1a']);
print(2,'-depsc',['example1b']);
print(3,'-depsc',['example1c']);
print(4,'-depsc',['example1d']);
print(5,'-depsc',['example1e']);
print(6,'-depsc',['example1f']);