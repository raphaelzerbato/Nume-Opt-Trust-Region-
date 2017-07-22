function varargout= trust_method_NM1(varargin)

%%% Fonction %%%
varin = varargin{1};
delta = varargin{2} 
dim = length(varin);
x = sym(zeros(1, dim));
for k=1:dim
    x(k) = sym(sprintf('x%d', k));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%FUNCTION TO CHANGE%%%%%%%%%%%%%%%%%%%%%%%%%%
F = x(1,1)^4+x(1,2)^2+9*(x(1,1)*x(1,2))^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


FX = gradient(F);

%%% calcule du gradient au point x_0 %%%

grad = zeros((dim),1);

for i=1:(dim)
    grad(i,1) = subs(FX(i,1), x, varin');
end

%%% calcule du Hessien %%%
H1 = hessian(F,x);
H = zeros(2,2);
for i=1:(dim);
    H(1,i) = subs(H1(1,i), x, varin');
    H(2,i)= subs(H1(2,i),x,varin');
end

%% construction de la fonction quadratic 
x_k = sym(zeros(1, dim));
for k=1:dim;
    x_k(k) = sym(sprintf('x%d', k));
end

quad_app =  f_2(varin)+ grad'*(x_k'-varin)+(1/2).*((x_k-varin')*H*(x_k'-varin));

%% Minimisation de Quad_app using Classical_Newton_M %%%
    
%%% Starting point of mini of Quad_app

[dkvect_CNTR] = CLASS_NEWTO_FUN_TR(varin,grad,H);
quad_app_image1 = double(subs(quad_app,x_k,varin'));
%%%% new point %%%
newcord_quad_app =varin+(delta*dkvect_CNTR);
quad_app_image2 =double(subs(quad_app,x_k,newcord_quad_app'));

%% Trust-region Method, verification du delta %%%
num= (f_2(varin)-f_2(newcord_quad_app));
den = (quad_app_image1-quad_app_image2);
R = num/den;
if num>0
    if  R<0.25 
        new_delta = delta/2;
    elseif R>0.75 
        new_delta= 1.5*delta;
    else
        new_delta = delta;
    end
else
    new_delta= delta/2;
    dkvect_CNTR = zeros(dim,1);
end
    %%% variable out

varargout{1} = new_delta;
varargout{2} = dkvect_CNTR;
varargout{3} = H;
  end

function varargout= CLASS_NEWTO_FUN_TR(varargin)

varin_CNTR = varargin{1};
grad_CNTR=varargin{2};
H_CNTR = varargin{3};
dimCNTR = length(varin_CNTR);

%%% Fonction %%%

mCNTR = dimCNTR; 
x_CNTR = sym(zeros(1, mCNTR));
for k=1:mCNTR
    x_CNTR(k) = sym(sprintf('x%d', k));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%FONCTION À CHANGER%%%%%%%%%%%%%%%%%%%%%%%%%%
F_CNTR = f_2(varin_CNTR)+ grad_CNTR'*(x_CNTR'-varin_CNTR)+(1/2).*((x_CNTR-varin_CNTR')*H_CNTR*(x_CNTR'-varin_CNTR));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FX_CNTR = gradient(F_CNTR);

%%% calcule du gradient au point x_0 %%%
grad_quad = zeros((dimCNTR),1);
grad1_quad = zeros((dimCNTR-1),1);
for i=1:(dimCNTR)
    grad1_quad(i,1) = subs(FX_CNTR(i,1), x_CNTR, varin_CNTR');
end

for i=1:(dimCNTR)
    grad_quad(i,1) = double(grad1_quad(i,1));
end

%%% calcule du Hessien %%%
H1_quad = hessian(F_CNTR,x_CNTR);
% H = double(H1);
H_quad = ones(2,2);
for i=1:(dimCNTR)
    H_quad(1,i) = subs(H1_quad(1,i), x_CNTR, varin_CNTR');
    H_quad(2,i)= subs(H1_quad(2,i),x_CNTR,varin_CNTR');
end

%%% steepest direction %%%
dkvect_CNTR = -(H_quad\grad_quad);

varargout{1} = dkvect_CNTR;
varargout{2} = H_quad;
end
