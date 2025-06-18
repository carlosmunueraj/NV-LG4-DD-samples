function [Ix, Iy, Iz, Id, Itotx, Itoty, Itotz] = operators(N)
px = [0 1; 1 0];
py = [0 -1i; 1i 0];
pz = [1 0; 0 -1];

for k = 1:N
    Ix(:, :, k) = kron(eye(2^(k-1)), kron(px, eye(2^(N-k))))/2;
    Iy(:, :, k) = kron(eye(2^(k-1)), kron(py, eye(2^(N-k))))/2;
    Iz(:, :, k) = kron(eye(2^(k-1)), kron(pz, eye(2^(N-k))))/2;
end

Id = eye(2^N);

Itotx = sum(Ix, 3);
Itoty = sum(Iy, 3);
Itotz = sum(Iz, 3);
end


