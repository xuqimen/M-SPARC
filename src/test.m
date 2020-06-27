function test()
% Mesh2ChebDegree(0.2)
% Mesh2ChebDegree(0.3)
% Mesh2ChebDegree(0.4)
% Mesh2ChebDegree(0.5)
% Mesh2ChebDegree(0.6)
% Mesh2ChebDegree(0.7)
% Mesh2ChebDegree(0.8)
% % 
% % 
% h = linspace(0.1,1.0);
% for i = 1:length(h)
% 	hi = h(i);
% 	npl(i) = Mesh2ChebDegree(hi);
% end
% for i = 1:length(h)
% 	hi = h(i);
% 	nplinvh(i) = Mesh2ChebDegreeInvh(hi);
% end
% 
% 
% figure; hold on; box on;
% plot(h,npl,'-','linewidth',2)
% plot(h,nplinvh,'--k','linewidth',2)
% %plot(h,11.2./h,'--r','linewidth',2)
% hold off;
% xlabel('h (Bohr)')
% ylabel('npl')
% %lg = legend('SPARC','$c/h$')
% lg = legend('SPARC','$P_m(1/h)$')
% set(lg,'interpreter','latex');
% 
% h = gcf;
% fname = 'npl'
% set(h,'Units','Inches');
% pos = get(h,'Position');
% set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[1.02*pos(3), 1.02*pos(4)])
% %saveas(gcf,fname,'epsc') 
% print(h,fname,'-dpdf','-r0')

% return


% test_ChebyshevSeries()


test_sqnpl();

end


function test_sqnpl()

Etotal_ref =-8.159022943;
force_ref = [   0.005167890961373   0.002128046516335   0.003809864488619
  -0.009367291979781   0.005801058180613   0.000190100551340
   0.009645473428618  -0.003131628251629   0.003146788524201
  -0.005446072410210  -0.004797476445320  -0.007146753564160];

% SQnplVec = [25:25:100,110:20:190,200:100:1000];
% SQnplVec = [110:20:190];
SQnplVec = 500;
E_err = zeros(size(SQnplVec));
F_err = zeros(size(SQnplVec));
nSCF  = zeros(size(SQnplVec));
ACTUAL_NPL = zeros(size(SQnplVec));

global mynpl
for i = 1:length(SQnplVec)
%for i = 1:2
	%fname = '../tests/MeshConvergence/Si8-ONCV-0.5'; 
	fname = '../mytest/Al4/sprc-calc';
	mynpl = SQnplVec(i)
	S = msparc(fname); 
	ACTUAL_NPL(i) = S.sq_npl;
	nSCF(i) = S.count_SCF;
	%E_err(i) = abs(S.Etotal - Etotal_ref)/S.n_atm;
	%F_err(i) = max(max(abs(S.force - force_ref)));
end

fprintf('npl #SCF  E_err  F_err \n');
for i = length(SQnplVec):-1:1
% 	fprintf('%4d  %.3e  %.3e\n',SQnplVec(i), E_err(i), F_err(i));
	%fprintf('%4d  %.3e  %.3e\n',ACTUAL_NPL(i), E_err(i), F_err(i));
	fprintf('%4d  %4d  %.3e  %.3e\n',ACTUAL_NPL(i),nSCF(i),E_err(i),F_err(i));
end

end


function test_ChebyshevSeries()
% x = linspace(-2,2)';
% y1 = chebyshevT(4,x); % using matlab inbuilt function (slow) 
% y2 = ChebyshevSum_matvec([0,0,0,0,1],diag(x),eye(length(x))); 
% plot(x,diag(y2),x,y1)

sigma = 0.003675;
%sigma = 0.001;
f = @(x) 1./(1+exp((x-0.24)/sigma));
%f = @(x) (x-1).^10 + 2*x;
%f = @(x) sin(x) + exp(x) + x.^7;
xmin = -0.5; xmax = 0.5;
% c = (xmax + xmin)/2;
% e = (xmax - xmin)/2;
% f_hat = @(y) f(e*y+c);
%f = @(x) chebyshevT(3,x);
npl = 100;
c = ChebyshevCoeff(npl, f, xmin, xmax);
disp(c)
close all;
figure; hold on; box on;
Ns = 100;
x = linspace(xmin,xmax,Ns);
tic; Y = ChebyshevSum_matvec(c,diag(x),eye(length(x)),xmin,xmax); toc; y = diag(Y); plot(x,y,'k-','linewidth',2)
plot(x,f(x),'r--','linewidth',2)
xlabel('$\lambda$','interpreter','latex')
ylabel('occ','interpreter','latex')
title(sprintf('npl = $%d$',npl),'interpreter','latex');
legend('Chebyshev series approx.','exact','location','southwest')
hold off;
h = gcf;
fname = 'smearing_fit'
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[1.02*pos(3), 1.02*pos(4)])
%saveas(gcf,fname,'epsc') 
print(h,fname,'-dpdf','-r0')

A = rand(Ns,Ns); A = (A+A')*0.5;
B = rand(Ns,Ns); B = B'*B;
tic
[C D] = eig(A,B);
toc

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function npl = Mesh2ChebDegree(h) 
	% the relation between h and npl is fit with a cubic polynomial
	% p(x) = p3 * x^3 + p2 * x^2 + p1 * x + p0.
	p3 = -700. / 3.;
	p2 = 1240. / 3.;
	p1 = -773. / 3.;
	p0 = 1078. / 15.;
	if (h > 0.7) 
		npl = 14;
	else 
		npl = ((p3 * h + p2) .* h + p1) .* h + p0;
	end
	npl = round(npl);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function npl = Mesh2ChebDegreeInvh(h) 
	% the relation between x = 1/h and npl is fit with a cubic polynomial
	% p(x) = p3 * x^3 + p2 * x^2 + p1 * x + p0.
	x = 1./h;
	p = [-112/1875, 242/375, 283/75, 112/15 ];
	
	if (h > 0.7) 
		npl = 14;
	else 
% 		npl = ((p3 * h + p2) .* h + p1) .* h + p0;
		npl = polyval(p,x);
	end
	npl = round(npl);
end