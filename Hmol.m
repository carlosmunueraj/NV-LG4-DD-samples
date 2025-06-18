function [Hcs, Hdd] = Hmol(nuc_coords, CS, euler, Ix, Iy, Iz, Bz, gh, mu0, hbar)

N = size(nuc_coords, 2);

% Rotations
Rz = @(x) [cos(x) -sin(x) 0; sin(x) cos(x) 0; 0 0 1];
Ry = @(x) [cos(x) 0 sin(x); 0 1 0; -sin(x) 0 cos(x)];
R = @(O) Rz(O(1))*Ry(O(2))*Rz(O(3));

% Rotate the nuclei to introduce variation
nuc_coords = R(euler)*nuc_coords % [0, pi/5, pi/3]

% Let us build the Hamiltonian
Hdd = 0;
Hcs = 0;

for k = 1:N
    Hcs = Hcs - (Bz * gh * CS(k) * 1e-6) * Iz(:, :, k); %Notice the (-) sign.   

    for l = (k+1):N
        rkl = nuc_coords(:, l) - nuc_coords(:, k);
        bkl = -(mu0 *hbar * gh^2)/(2 * norm(rkl)^3); % Dipole-dipole constant.
        
        rz = rkl(3)/norm(rkl);
        dkl = bkl * (3*rz^2 - 1)/2; % Secular dipole-dipole constant.
        
        Hdd = Hdd + dkl * (3*Iz(:, :, k)*Iz(:, :, l)-(Ix(:, :, k)*Ix(:, :, l)+...
              Iy(:, :, k)*Iy(:, :, l)+Iz(:, :, k)*Iz(:, :, l))); % Dipole-dipole hamiltonian.
    end
end

end

