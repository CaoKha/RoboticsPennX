boundary = [0; 2; 3; 1; 1; 5; 2; 3];
ps = [2];
ts = [1; 3];
% syms t0 t1 t21 t32 real
t = sym('t',[size(ts,1)+1 1]);
% position constraint in matrix form
nb_row_pos = 2*(size(ps,1)+1);
nb_col_pos = 8*(size(ps,1)+1);
D_pos = sym(zeros(nb_row_pos,nb_col_pos));
k = 1;
p = 1;
for i = 1:2:nb_row_pos
        m = 0;
        for j = k:k+7
            D_pos(i,j) = t(1)^m;
            D_pos(i+1,j) = t(p+1)^m;
            m = m + 1;
        end
    k = k + 8;
    p = p + 1;
end


nb_row_deriv = 3*(size(ts,1)+1)+3*size(ps,1);
nb_col_deriv = nb_col_pos;
D_deriv = sym(zeros(nb_row_deriv,nb_col_deriv));
k = 1;
p = 1;
q = 0;
for i = 1:3:nb_row_deriv
    if i == 1 
        m = 0;
        for j = k:k+7
            D_deriv(i,j) = diff(t(1)^m);
            D_deriv(i+1,j) = diff(D_deriv(i,j));
            D_deriv(i+2,j) = diff(D_deriv(i+1,j));
            m = m+1;
        end
    elseif i == (4 + q*6) && i ~= (nb_row_deriv-2)
        m = 0;
        for j = k:k+7
            D_deriv(i,j) = diff(t(p+1)^m);
            D_deriv(i+1,j) = diff(D_deriv(i,j));
            D_deriv(i+2,j) = diff(D_deriv(i+1,j));
            D_deriv(i+3,j) = diff(D_deriv(i+2,j));
            D_deriv(i+4,j) = diff(D_deriv(i+3,j));
            D_deriv(i+5,j) = diff(D_deriv(i+4,j));
            m = m+1;
        end
        m = 0;
        for j = k+8:k+15
            D_deriv(i,j) = -diff(t(1)^m);
            D_deriv(i+1,j) = diff(D_deriv(i,j));
            D_deriv(i+2,j) = diff(D_deriv(i+1,j));
            D_deriv(i+3,j) = diff(D_deriv(i+2,j));
            D_deriv(i+4,j) = diff(D_deriv(i+3,j));
            D_deriv(i+5,j) = diff(D_deriv(i+4,j));
            m = m + 1;
        end
        k = k + 8;
        p = p + 1;
        q = q + 1;
    elseif i == (nb_row_deriv-2)
        m = 0;
        for j = k:k+7
            D_deriv(i,j) = diff(t(end)^m);
            D_deriv(i+1,j) = diff(D_deriv(i,j));
            D_deriv(i+2,j) = diff(D_deriv(i+1,j));
            m = m+1;
        end
        break;
    end
end

q_pos = zeros(nb_row_pos,1);
h = 1;
for i = 1:2:nb_row_pos
    if i == 1
        q_pos(i) = boundary(1);
        q_pos(i+1) = ps(h);
    elseif i == (nb_row_pos-1)
        q_pos(i) = ps(h);
        q_pos(i+1) = boundary(5);
    else
        q_pos(i) = ps(h);
        h = h +1;
        q_pos(i+1) = ps(h);
    end
end

q_deriv = zeros(nb_row_deriv,1);
for i = 1:nb_row_deriv
    if i == 1
        q_deriv(i:i+2) = boundary(2:4);
    elseif i == (nb_row_deriv - 2)
        q_deriv(i:i+2) = boundary(6:8);
        break;
    elseif i > 3 
        q_deriv(i) = 0;
    end
end
big_D = [D_pos; D_deriv];
big_q = [q_pos; q_deriv];

a = zeros(size(ts,1)+1,1);
for i = 1:size(t,1)
    if i == 1
        a(i) = 0;
    elseif i == 2
        a(i) = ts(1);
    else 
        a(i) = ts(i-1) - ts(i-2);
    end
end
big_D = subs(big_D, t(1:end), a(1:end));
C_col = big_D\big_q;
num_col = size(C_col,1)/8;
C = zeros(8,num_col);
k = 1;
for i = 1:num_col
    C(1:8,i) = C_col(k:k+7);
    k = k + 8;
end
plot_traj(C,ts); 


