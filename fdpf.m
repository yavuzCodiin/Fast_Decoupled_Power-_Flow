tic;                                                  % starts a stopwatch timer.
clear;                                                % clears all variables from the workspace.
clc;                                                  % clears the command window.
load Case9ybus;                                       % loads data from the file 'Case9ybus.mat' into the workspace.
rbusdata = importdata('case9bus.xlsx');               % imports data from the Excel file 'case9bus.xlsx' and stores it in the variable rbusdata
busno = 9;                                            % initializes the variable busno with a value of 9.
 
 
branchdata = importdata('case9branch.xlsx');          % imports data from the Excel file 'case9branch.xlsx' and stores it in the variable branchdata
 
 
fbus = branchdata(:,1);                               % extracts the first column of branchdata and stores it in the variable fbus
fbusno = length(fbus);                                % calculates the number of elements in fbus and stores it in fbusno
branchno = fbusno-1;                                  % calculates the number of branches (branchno) by subtracting 1 from fbusno
%These lines extract specific columns from branchdata and assign them to variables fbus, tbus, r, x, b, a. It also calculates z and y based on extracted data.
fbus = branchdata(1:branchno,1);                       
tbus = branchdata(1:branchno,2);                      
r = branchdata(1:branchno,7);                          
x = branchdata(1:branchno,8);                         
b = 1i*branchdata(1:branchno,9);                      
a = branchdata(1:branchno,10);
z = r + x*1i; 
y = 1./z; 
%These lines initialize Ybus as a zero matrix and then populate it based on the values of fbus, tbus, y, a, and b calculated earlier.
Ybus = zeros(busno,busno); 
 for k = 1:branchno 
 if a(k)== 0 
 a(k)= 1 ; 
 end
 end
 
for i = 1:busno 
 for j = 1:branchno 
 if fbus(j) == i 
 YBus(i,i) = YBus(i,i) + y(j)/(a(j)^2) + b(j); 
 elseif tbus(j) == i 
 YBus(i,i) = YBus(i,i) + y(j) + b(j); 
 end
 end
end
%These lines import additional data from 'case9bus.xlsx' and extract relevant columns.
busdata = importdata('case9bus.xlsx');                % BusType categorizes buses into different types (slack, generator, load).          
BusType = busdata(:,3);                               % the third column of the busdata matrix to determine the Type of the Bus.
slackbus = find(BusType==3);                          % slackbus stores the indices of slack buses.
PVbus = find(BusType==2 | BusType==3);                % PVbus stores indices of PV (generator) buses.

mgen = length(PVbus); 
PQbus = find(BusType==0 | BusType==1);                % PQbus stores indices of PQ (load) buses.
loadnumber = length(PQbus); 
%Pd, Qd, Pg, and Qg store real and reactive power values for loads and generators, converted to per unit (p.u.) values.
Pd = busdata(:,6);  
Pd = Pd./100; 
Qd = busdata(:,7); 
Qd = Qd./100; 
Pg = busdata(:,8);  
Pg(1:1) = 0; 
Pg = Pg./100; 
Qg = busdata(:,9); 
Qg(1:1) = 0; 
Qg = Qg./100;
%Power Flow 
%P and Q are calculated as the differences between generator and load values.
P = Pg - Pd; % Pi = Pgi -
Pdi
Q = Qg - Qd; % Qi = Qgi -
Qdi
Ppf=P; 
Qpf=Q; 
V = busdata(:,11); 
%These lines prepare matrices for calculation (G, B, Brp, Brq, Brpinv, Brqinv) needed for the power flow calculations.
for u = 1:busno 
 if V(u)== 0 % If it is not zero it measns it is from generator buses, if it is zero 
 V(u)= 1; % this makes them 1 for the flat start.
 end
end
 angle = zeros(busno,1); 
 % Ybus = G + j*B
G = real(YBus); 
B = imag(YBus); 
 % Br,p 
 
Brp = B;
Brp(:,1) = [];
Brp(1,:) = [];
Brpinv = inv(Brp); %inverse of Br,p. 
% Br,q
 
Brq = B;
for i=1:mgen
 Brq(:,PVbus(i)-i+1) = [];
 Brq(PVbus(i)-i+1,:) = []; 
end
Brqinv = inv(Brq); %inverse of Br,q. 

%This section calculates power mismatches and iteratively adjusts voltage angles and magnitudes until the mismatches are 
%below the specified threshold (0.001) or the maximum iterations (10) are reached. The rest of the code deals with 
%post-processing and output formatting. It calculates line currents, line power flows, and prints out the results.

P_V = (Pg-Pd)./abs(V); % Pi = Pgi - Pdi
Q_V = (Qg-Qd)./abs(V); % Qi = Qgi - Qdi
Psp = P_V;
Qsp = Q_V;
Epsilon=1;
Iter=0;
while (Epsilon > 0.001 && Iter<10) 
 
 P_V = zeros(busno,1);
 Q_V = zeros(busno,1);
 
 % P and Q calculations
 for i = 1:busno
 for k = 1:busno
 P_V(i) = P_V(i) + V(k)*(G(i,k)*cos(angle(i)-
angle(k)) + B(i,k)*sin(angle(i)-angle(k))); % Real Power Calculation
 Q_V(i) = Q_V(i) + V(k)*(G(i,k)*sin(angle(i)-
angle(k)) - B(i,k)*cos(angle(i)-angle(k))); % Reactive Power Calculation
 end
 end
 
 
 fPn = P_V-Psp; 
 fQn = Q_V-Qsp; 
 fP = fPn(2:busno); % calculates the unknown real power values matrix
 k = 1;
 fQ = zeros(loadnumber,1);
 for i = 1:busno
 if BusType(i) == 0
 fQ(k,1) = fQn(i); % calculates the unknown reactive power values matrix
 k = k+1;
 end
 end
 
 
 f = [fP; fQ]; % It calculates the f(x) matrix 
 
 dangle = Brpinv * fP; 
 dvolt = Brqinv .* fQn;
 angle(2:busno)= angle(2:busno)+ dangle;
 Iter = Iter + 1; 
 Epsilon = max(abs(f)); 
 
end % When Epsilon is less than 0,001 iteration finishes.
 
Vcomplex = diag(V)*cos(angle) + 1j*diag(V)*sin(angle); 
 % the Line Currents:
Iij = zeros(busno,busno); 
Sij = zeros(busno,busno); 
for m = 1:branchno
 p = fbus(m); 
 q = tbus(m);
 
end
Iij = sparse(Iij); 
 % Line Power Flows:
for m = 1:busno
 for n = 1:busno
 if m ~= n
 Sij(m,n) = Vcomplex(m)*conj(Iij(m,n))*100;
 end
 end
end
Sij = sparse(Sij);
Pij = real(Sij);
Qij = imag(Sij);
angle = (angle*180)/pi; 
 % Program Outpus:
fprintf('\n');
fprintf('*** Outputs***:\n');
fprintf('\n');
fprintf('\n');
display (YBus) 
fprintf('\n*Bus Data Information:\n')
fprintf(' \tVOLTAGE\t\tANGLE(deg)\tPin (MW) \tQin 
(MVar)\n');
fprintf('Bus %d: %8f \t%f\t%f\t%f\n',[(1:busno)' V angle 
P*100 Q*100]'); 
fprintf('\n');
fprintf('Number of iterations: %d \n',Iter);
fprintf('\n');
disp('Number of Buses:')
fprintf('%d\n', busno);
toc;
