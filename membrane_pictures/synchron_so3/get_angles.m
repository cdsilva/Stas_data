function A = get_angles(R)

A = zeros(1, 3);

A(1) = atan2(R(2,3),R(1,3));
A(2) = acos(R(3,3));
A(3) = atan2(R(3,2), -R(3,1));

A = real(A);

