
q = [2,0,9,2,0,3,7,5,1,2,0,7,6,9,1,2,6,2]; %ID1_ID2
ro=1;
M=18;%matrix size
h = [1/5,1/2,2,5,10];
q = q';


%initializing
arr_det = zeros(5,1);
arr_error_compute = zeros(5,1);


for i = 1:5
%fillin the matrix    
A_3 = fill_A(zeros(M),h(i));
v_tag_3 = A_3*q; % calculating matrix multiplaction

arr_det(i) = abs(det(A_3));% calculate the Matrix determinant


kappa_A_3 = cond(A_3,"inf"); %finding condition number
disp ("number of condition number: "+ kappa_A_3 ); %print

%making pseudo invers matrix
psi = pinv(A_3); %(ATA)^-1*AT

close_q = psi * v_tag_3;
minimum_last_squre = norm(v_tag_3 - A_3 * close_q,2);
disp("the minimum lasp squre is: " + minimum_last_squre )
arr_error_compute(i) = abs(norm(q - close_q,2))./norm(q,2);

end


%-------------ploting------------
figure('name',"question 3")
plot_3 = loglog(h,arr_det,h,arr_error_compute);
legend("determinant","error |q-q'|");
plot_3(1).LineWidth = 3; %Change width of the line in the graph
plot_3(2).LineWidth = 3;





%-----------------------functions -----------------------
function A = fill_A(A,h_co)
M = size(A,1);%numbers of electrostati charges - q 
ro =1;
h = h_co* pi.* ro./ M;
r=0;
for m = 1:M
    for n = 1:M
        r = sqrt((h+ro*sin(((m*pi)/M))-ro*sin(((n*pi)/M))).^2+(ro*cos((m*pi)/M)-ro*cos((n*pi)/M)).^2);
        formula = 4*pi.*r;
        A(m, n) = 1./formula;
    end
end
end
