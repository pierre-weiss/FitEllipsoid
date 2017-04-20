function q=Project_on_B(q0)
Q0=[[q0(1),q0(4)/sqrt(2),q0(5)/sqrt(2)];[q0(4)/sqrt(2),q0(2),q0(6)/sqrt(2)];...
    [q0(5)/sqrt(2),q0(6)/sqrt(2),q0(3)]];
[U,S0]=eig(Q0);
s0=diag(S0);
s=projsplx(s0);
S=diag(s);
Q=U*S*U';
q=[Q(1,1);Q(2,2);Q(3,3);sqrt(2)*Q(2,1);sqrt(2)*Q(3,1);sqrt(2)*Q(3,2);q0(7:end)];
end