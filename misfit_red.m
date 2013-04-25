function [f,g] = misfit_red(mt,Q,D,model)
% reduced objective:
%
%  .5*||PA(m)^{-1}q - d||^2.
%
% use:
%   [f,g] = misfit_red(mt,Q,D,model)
%
% input:
%   mt - model [s^2/m^2]
%   Q  - sources
%   D  - data
%   model.h - [dz,dx] gridspacing in z and x direction [m]
%   model.n - [nz,nx] number of gridpoints in z and x direction
%   model.f - frequencies
%   model.Is - source location indices
%   model.Ir - receiver location indices
%
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

% get sampling operators
Ps = getP(model.n,model.Is{:});
Pr = getP(model.n,model.Ir{:});

% initialize misfit and gradient
f = 0;
g = zeros(prod(model.n),1);

% loop over frequencies 
for k = 1:length(model.f)
	% get Helmholtz operator
	At = getA(model.f(k),mt(:),model.h,model.n);

	% solve forward wave-equation and sample data
	Ut = At\Ps'*Q;
	Dt = Pr*Ut;

	% solve adjoint wave-equation
	Vt = At'\(Pr'*(Dt - D));

	% compute misfit and gradient
	f = f + .5*norm(Dt(:) - D(:)).^2;
	g = g -(2*pi*model.f(k))^2*real(sum(conj(Ut).*Vt,2));
end