% the Anitescu-Potra model with Newtonian restitution
%
% @article{Anitescu:1997,
%	Author = {Miha Anitescu and Florian A. Potra},
%	Journal = {Nonlinear Dynamics},
%	Pages = {231--247},
%	Title = {Formulating Dynamic Multi-Rigid-Body Contact Problems with Friction as Solvable Linear Complementarity Problems},
%	Volume = {14},
%	Year = {1997}}

function [v_plus, z] = APNewton(M, n, s, v, ha, mu, e)
    
   global quadProgOpts
   
   E = [1;1];
   D = [-s, s];
   
   % Define LCP:
   M_hat = [n'*(M \ n),  n'*(M \ D),  0;
            D'*(M \ n),  D'*(M \ D),  E;
                    mu,         -E',  0];

   % Nima: note that you are making a modeling choice here by not
   %       using the old (the not integrated) velocity 
   q_hat = [n'*(v + ha) + n'*e*v; D'*(v + ha); 0];
   
   withLemke = 0;
   
   if withLemke
       z = LCP(M_hat,q_hat);
   else
%        options=optimoptions('quadprog',...
%                 'Algorithm','interior-point-convex','Display','off'); 
       A = [-M_hat; -eye(size(M_hat,1))];
       B = [ q_hat;zeros(size(M_hat,1),1)];
       M_hat = (M_hat+M_hat')/2;
%        z = quadprog(2*M_hat,q_hat,A,B,[],[],[],[],[],options);
       z = quadprog(2*M_hat,q_hat,A,B,[],[],[],[],[],quadProgOpts);
   end
   
   %         if abs(z(1))<0.05
   %             e = 0;
   % %             vlplus1 = zeros(3,1)+1e-3;
   %         end
   
   v_plus = [M \ n, M \ D]*[z(1:3)]+(v + ha);

