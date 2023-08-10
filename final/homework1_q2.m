%% 
%{
read me!: use section running to operate the two scrip: first:
 section 1:initialize general variables

 section 2: gauss-seddel script and graph plot;

 section 3: jacobi script and graph plot;
%}

%------------------------------question 2---------------------------

q = [2,0,9,2,0,3,7,5,1,2,0,7,6,9,1,2,6,2]';
ro=1;
M=18;
A_2 =zeros(M);
A=zeros(18);
iteration_num = 0;
%% 

for i=[5,2,1]
A_2 = fill_A(A_2,1,i);  
iteration_num = 1;

 
k_axis = [];
iter_axis = [];
real_err =[];

v_tag_2 = A_2*q;


%---------------------------------gausse seidell---------------------------

  
L = tril(A_2,-1); %lower part of marix
D = diag(diag(A_2));    %diagonal part of marix
U = triu(A_2,1); %upper part of marix
Q = L + D;       % Q matrix
Q_inv = inv(Q);  %Q^-1

C = Q_inv * v_tag_2;
q_k =  C;   %initializ q_k when given start q_0 = 0
G = -Q_inv*U;

G_norm = norm(G,"inf");

while(abs(norm(q-q_k,"inf")./norm(q,"inf")) > 10^-3)    
q_k_min_1 = q_k; 

%real relative error caculate
real_err = [real_err , norm(q_k-q,"inf")./norm(q,"inf")];

 % iteration step
q_k = G * q_k_min_1 + C;
iter_axis =[iter_axis, iteration_num] ; 

%calculate relative errore between q_k-q_k-1
k_axis =[k_axis, norm(q_k-q_k_min_1,"inf")./norm(q_k_min_1,"inf")] ;

 %increas iteration step
iteration_num = iteration_num + 1;

%self bound for converjion not matematicali form
if iteration_num>10000
    disp("not converge")
end
end

xxx = norm(q-q_k,"inf")./norm(q,"inf");


%graph ploting
figure('name', "gause seidell "+ i);
graph = semilogy(iter_axis,k_axis,iter_axis,real_err);
graph(1).LineWidth = 1; %Change width of the line in the graph
graph(2).LineWidth = 1;
legend("q(k)-Q(k-1) rel dis","real error");
title("gauusi-seidel:relative error iteration and real error in func of iteration num: M*"+ i );
xlabel("iteration number");
ylabel("error");


end 
%% 
%----------------------------jacobi------------------------------------

%-------------------------question 2.c-------------------------


%making A_2 matrix and fill it

k_axis = [];
iter_axis = [];
real_err =[];
A_2 = zeros(M);
A_2 = fill_A(A_2,2,5);
v_tag_2 = A_2*q; 
iteration_num = 1;

%checking if dominant diagonal 
is_dominant(A_2);


%making inverse Q matrix 
Q = diag(diag(A_2));
invQ = zeros(18);
for i = 1:M
    invQ(i,i) = 1./ Q(i,i);
end

%making iteration formula
I = eye(M);
G = I -invQ*A_2;
G_norm = norm(G,"inf");    
C = invQ*v_tag_2;  
q_k = C;

%cheking if it converge
if G_norm <1
    disp("there is conversion!");
end


%step iteration
while(abs(norm(q-q_k,"inf")/norm(q,"inf")) > 10^-3) %if relativ < tollerance   
q_k_min_1 = q_k;
real_err = [real_err , norm(q_k-q,"inf")./norm(q,"inf")];

%q_k = G*q_k-1 +Q^-1*v
q_k = G * q_k_min_1 + C; 

%inser data of iteration
iter_axis =[iter_axis, iteration_num] ; 
k_axis =[k_axis, norm(q_k-q_k_min_1,"inf")./norm(q_k_min_1,"inf")] ;
iteration_num = iteration_num + 1;
end




%grapg plotting
figure('name',"yakobi");graph_2 = semilogy(iter_axis,k_axis,iter_axis,real_err);
graph_2(1).LineWidth = 3; %Change width of the line in the graph
graph_2(2).LineWidth = 3;
legend("q(k)-Q(k-1) rel dis","real error");
title("jacobi-relative error iteration and real error in func of iteration num: M*"+ 5 );
xlabel("iteration number");
ylabel("error");


%fill_A_matrix
function A = fill_A(A,x,i) %x-number of degree of r. i-coefficient
M = size(A,1);%numbers of electrostati charges - q 
ro =1;
h = pi.* ro./ (M.*i);
r=0;
for m = 1:M
    for n = 1:M
        r = sqrt((h+ro*sin(((m*pi)/M))-ro*sin(((n*pi)/M))).^2+(ro*cos((m*pi)/M)-ro*cos((n*pi)/M)).^2);
        formula = 4*pi.*r^x;
        A(m, n) = 1./formula;
    end
end
end

%--------------ckecking dominant 
 function flag = is_dominant(A)
 d = diag(A);
 flag =false;
 M =size(A,1);
 for i= 1:M
     sum = 0;
     for j = 1:M
         if i ~= j
         sum = sum + abs(A(i,j));
         end
     
     end
     if sum >= abs(d(i))
         disp("not diminant diagonal in qution");
         flag = true;
         break;
     end
     sum = 0;
 end
 end


