% electrostati charges - q
M =18; %numbers of electrostati charges - q 
q = [2,0,9,2,0,3,7,5,1,2,0,7,6,9,1,2,6,2];


q_list = zeros(6,18);
%r = rho(m)
r=1;
%coefficient    
coeff = [1,2,5,10,20,50];

%initialize arr errirs each letter indicates the number of question part
coeff_plot = coeff.*((pi * r) ./ M);

arr_err_mat_A = zeros(6);
arr_err_B = zeros(6,1); % error |q-q_tag|
arr_err_C = zeros(6,1);% error |q-q_tag| Aq=v+delta*v
arr_err_D = zeros(6,1);% error |q-q_tag| (A+delta*A)*q=v   delta <<1
arr_kappa = zeros(6,1); %condiotn number of A
arr_q_norm = zeros(6,1); %q_norm |q|
arr_v_tag_norm =zeros(6,1);% v norm |v|
arr_A_norm = zeros(6,1);% A norm |A|    
for c = 1:6
%height of center of two paraller circle arc
co = coeff(c); 
h = co*pi*r/M; 


%make matrix of zeros size M
A = zeros(M);



% Populate the matrix with values
for m = 1:M
    for n = 1:M
        A(m, n) = 1./(4*pi*sqrt((h+r*sin(((m*pi)/M))-r*sin(((n*pi)/M))).^2+(r*cos((m*pi)/M)-r*cos((n*pi)/M)).^2));
    end
end


% calculat A * q
v_tag = A*q';


%lu calculate
[L, U, P] = lu(A);

%Condition number calculate
arr_kappa(c) = cond(A,"inf");

%norm calculate
   
arr_q_norm(c) = norm(q,2);

arr_v_tag_norm(c) = norm(v_tag,2);

arr_A_norm(c) = norm(A,"fro");






      







%---------------------------------------------------question B---------


%calculating Aq' = v'
y = ly(L,P*v_tag,M);
q_tagi = Ux(U,y,M);
q_list(c,:) = q_tagi;
%relative error
delta = norm(q_tagi-q',2);
error = delta ./ arr_q_norm(c);
arr_err_B(c)=error;

%---------------------------------------------------question c---------
delta_v =10.^-3.*arr_v_tag_norm(c);
v_new = v_tag+delta_v;

%calculating Aq' = v'+delta v_m
y_new = ly(L,P*v_new,M);
q_tag_new = Ux(U,y_new,M);
delta_new = norm(q_tag_new-q',2)./norm(q,2);
arr_err_C(c) = norm(q_tag_new-q',2)./norm(q,2);

%---------------------------------------------------question D---------
delta_A =10^-3.*arr_A_norm(c);
A_new = A + delta_A;

%calculating A_new*q' = v'
[L_new, U_new, P_new] = lu(A_new);
y_new_D = ly(L_new,P_new*v_tag,M);
q_tag_new_D = Ux(U_new,y_new_D,M);
delta_new_D = norm(q_tag_new_D-q',2)./norm(q,2);
arr_err_D(c) = norm(q_tag_new_D-q',2)./norm(q,2);


end
%-------------------------------ploting--------------------------------
figure('name',"question 1")
plt=loglog(coeff_plot,arr_kappa,"-*",coeff_plot,arr_err_B,"-*",coeff_plot,arr_err_C,"-*",coeff_plot,arr_err_D);
plt(1).LineWidth = 1; %Change width of the line in the graph
plt(2).LineWidth = 1;
plt(3).LineWidth = 1;
plt(4).LineWidth = 1;

legend("cond num","RL B","RL C","RL D");
xlabel("distance");
ylabel("Realtive error and condition number");
title("Realtive error and condition number");
grid on;

%{
%------------------------------question 2---------------------------

h = pi.* r./ (5.* M);
A_2 = zeros(M);
for m = 1:M
    for n = 1:M
        A_2(m, n) = 1./(4*pi*sqrt((h+r*sin(((m*pi)/M))-r*sin(((n*pi)/M))).^2+(r*cos((m*pi)/M)-r*cos((n*pi)/M)).^2));
    end
end

v_tag_2 = A_2*q';
q_tag = zeros(18,1);
[L,U] = lu(A_2);
D = diag(A_2);
iteration_num = 0;
k_axis = [];
iter_axis = [];
real_err =[];
while(norm(q-q_tag_2,2) < 10^-3)    
Q = L + D;
Q_inv = inv(Q);
q_last = q_tag_2;
C = Q_inv\v_tag_2;
real_err = [real_err,norm(q_tag_2-q,"inf")./norm(q,"inf")];
q_tag_2 = -Q_inv*U*q_tag_2 + C;
iteration_num = iteration_num + 1;
iter_axis =[iter_axis, iteration_num + 1] ;
k_axis =[k_axis, norm(q_tag_2-q_last,"inf")./norm(q_last,"inf")] ;
end

semilogy(iter_axis,k_axis,real_err,k_axis);



%}
%------------------------------functions---------------------------

%calculate Ly = b
function y = ly(L,b,M)
    y = zeros(M,1);
    y(1) = b(1)/L(1,1);
    for i=2:M
        accumulate = 0;
        for j=1:i-1
            accumulate = accumulate + L(i,j).*y(j);
        end

        y(i) =(b(i) - accumulate)./ L(i,i);
    end
end





%calculate Ux = y

function x = Ux(U,y,M)
    x = zeros(M,1);
    x(M) = y(M) ./ U(M,M);
    for i= M-1:-1:1
        accumulate = 0;
        for j=i+1:M
            accumulate = accumulate + U(i,j).*x(j);
        end
        x(i) = (y(i) - accumulate)./ U(i,i);
    end
end

%{
%pivoteing function - return a permutation matrix of eye
function p = pivoting(A,idx)
    i = eye(18);
    [row,~] =find(max(A(:, idx)));
    p(row,:) = i(0, : );
    p(0, :) =  i(row, :);
end
%}



%{

A_iter = A;
l = eye(18);
u = eye
for i = 1:M-1
    p = pivoting(A_iter,i);
    A_iter = p*A_iter;
    b = p*b;
    l = eye(18);
    for j = i+1:M-1
        l(j,i) = -A_iter(j,i)/A_iter(i,i);
    end

    A_iter = l*A_iter;
end

   %}








   




