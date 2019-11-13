syms c1 c2 c3 c4 c5 c6 s1 s2 s3 s4 s5 s6
a = 13;
b = 2.5;
c = 8;
d = 2.5;
e = 8;
f = 2.5;
T01 = [c1 0 s1 0;
    s1 0 -c1 0;
    0 1 0 a;
    0 0 0 1];
T12 = [c2 -s2 0 c*c2;
    s2 c2 0 c*s2;
    0 0 1 -b;
    0 0 0 1];
T23 = [c3 0 -s3 0;
    s3 0 c3 0;
    0 -1 0 -d;
    0 0 0 1];
T34 = [c4 0 s4 0;
    s4 0 -c4 0;
    0 1 0 e;
    0 0 0 1];
T45 = [c5 0 -s5 0;
    s5 0 c5 0;
    0 -1 0 0;
    0 0 0 1];
T56 = [c6 -s6 0 0;
    s6 c6 0 0;
    0 0 1 f;
    0 0 0 1];
T02 = T01*T12;
T03 = T02*T23;
T04 = T03*T34;
T05 = T04*T45;
T06 = T05*T56;
T36 = T34*T45*T56;