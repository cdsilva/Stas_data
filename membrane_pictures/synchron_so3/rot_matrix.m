function R = rot_matrix(A)
R = [cos(A(1))*cos(A(2))*cos(A(3))-sin(A(1))*sin(A(3))  -cos(A(3))*sin(A(1))-cos(A(1))*cos(A(2))*sin(A(3)) cos(A(1))*sin(A(2)) ;
    cos(A(1))*sin(A(3))+cos(A(2))*cos(A(3))*sin(A(1)) cos(A(1))*cos(A(3))-cos(A(2))*sin(A(1))*sin(A(3)) sin(A(1))*sin(A(2));
    -cos(A(3))*sin(A(2)) sin(A(2))*sin(A(3)) cos(A(2))];