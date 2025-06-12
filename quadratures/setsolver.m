function [s,A,c,b,eiginvA,V,invV]=setsolver(solver)

%Set Runge-Kutta solver
if strcmp(solver,'LobattoIIIC2')
    s = 2;
    c = [0  1]';
    b = [1/2  1/2]';
    A = [1/2,  -1/2;
        1/2,   1/2];
    [V,D] = eig(A);
    invV = inv(V);
    eiginvA=1./diag(D);
    
elseif strcmp(solver,'LobattoIIIC4')
    s = 3;
    c = [0      1/2    1]';
    b = [1/6    2/3    1/6]';
    A = [1/6,  -1/3,   1/6;
        1/6,   5/12, -1/12;
        1/6,   2/3,   1/6];
    [V,D]=eig(A);
    invV = inv(V);
    eiginvA=1./diag(D);
    
elseif strcmp(solver,'LobattoIIIC6')
    sq5 = sqrt(5); s = 4;
    c = [0      (5-sq5)/10          (5+sq5)/10          1]';
    b = [1/12    5/12                5/12            1/12]';
    A = [1/12,   -sq5/12,            sq5/12,        -1/12;
        1/12,   1/4,                (10-7*sq5)/60,  sq5/60;
        1/12,   (10+7*sq5)/60,      1/4,           -sq5/60;
        1/12,   5/12,               5/12,           1/12];
    [V,D]=eig(A);
    invV = inv(V);
    eiginvA=1./diag(D);
    
elseif strcmp(solver,'BE')||strcmp(solver,'RadauIIA1')||solver==1
    s = 1; c=1; b=1; A=1;
    V=1; invV = 1;
    eiginvA=1;
    
elseif strcmp(solver,'RadauIIA3')||solver==2
    s = 2; c = [1/3 1]'; b = [3/4 1/4]';
    A  = [5/12 -1/12;
        3/4  1/4];
    [V,D]=eig(A);
    invV = inv(V);
    eiginvA=1./diag(D);
    
elseif strcmp(solver,'RadauIIA5')||solver==3
    sq6 = sqrt(6); s = 3;
    c = [(4-sq6)/10 (4+sq6)/10 1]';
    b = [(16-sq6)/36 (16+sq6)/36 1/9 ]';
    A = [(88-7*sq6)/360,        (296-169*sq6)/1800,    (-2+3*sq6)/225 ;
        (296+169*sq6)/1800,     (88+7*sq6)/360,        (-2-3*sq6)/225 ;
        (16-sq6)/36,            (16+sq6)/36,           1/9             ];
    
    [V,D]=eig(A);
    invV = inv(V);
    eiginvA=1./diag(D);
    
end
