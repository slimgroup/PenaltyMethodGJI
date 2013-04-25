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
n  = [51 101];
h  = [20 20];
z  = [0:n(1)-1]'*h(1);
x  = [0:n(2)-1]*h(2);
[zz,xx] = ndgrid(z,x);

v0 = 2000*ones(n);
dv = zeros(n);
dv = dv - 100*(abs(zz-200)<=10);
dv = dv + 200*((zz-400).^2+ (xx-1000).^2<=1000);
dv = dv + 200*(abs(zz+.25*xx-1000)<=10);

% initialized images
I1 = zeros(prod(n),1);
I2 = zeros(prod(n),1);

% loop over frequencies
for f = [2:10]
    % define operators
    Ps = getP(n,2,[1:101]);
    Pr = getP(n,2,[1:101]);
    A  = getA(f,1./(v0(:) + dv(:)).^2,h,n);
    A0 = getA(f,1./(v0(:)).^2,h,n);
    Q  = speye(101);

    % make data
    U  = A\(Ps'*Q);
    D  = Pr*U;

    % RTM
    U1 = A0\(Ps'*Q);
    V1 = A0'\(Pr'*(D - Pr*U1));
    I1 = I1 + real(sum(conj(U1).*V1,2));

    % penalty method
    lambda = 1;
    U1 = [lambda*A0;Pr]\[lambda*Ps'*Q;D];
    I2 = I2 + real(sum(conj(U1).*(A0*U1-Ps'*Q),2))./((2*pi*f)^2*sum(abs(U1).^2,2));
end

% plot results
figure;imagesc(x,z,dv,[-1 1]*2e2);colormap(gray);axis equal tight
xlabel('x [m]','fontsize',20);ylabel('z [m]','fontsize',20);
set(gca,'fontsize',20);

figure;imagesc(x,z,reshape(I1/norm(I1,'inf'),n),[-1 1]);colormap(gray);axis equal tight
xlabel('x [m]','fontsize',20);ylabel('z [m]','fontsize',20);
set(gca,'fontsize',20);

figure;imagesc(x,z,reshape(I2,n)/norm(I2,'inf'),[-1 1]);colormap(gray);axis equal tight
xlabel('x [m]','fontsize',20);ylabel('z [m]','fontsize',20);
set(gca,'fontsize',20);

% save figures
print(1,'-depsc',['example4a']);
print(2,'-depsc',['example4b']);
print(3,'-depsc',['example4c']);
