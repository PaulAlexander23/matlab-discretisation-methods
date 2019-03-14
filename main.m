
n = 2^8;
x = linspace(1/n,1,n)';

y = cos(2*pi*x);

problem = "fTest";
problemDeg = [1,2];

method = "diff_ps";

if method == "diff_fd"
    init_fd(x, 2);
end

tic
f = feval(problem,x,y,method);
tc = toc;

if method == "diff_fd"
    del_fd();
end

fexact = -sin(2*pi*x);

fe = f - fexact;

err = norm(fe,'inf');

fprintf('Problem: %s, Method: %s, Error: %g, Time: %gs\n',problem,method,err,tc);

plot(x,f,x,fexact);

%%

n = [2^7,2^8];
L = [1,2];
x = cell(2,1);
for ni = 1:length(n)
    x{ni} = linspace(1/n(ni),1,n(ni))'*L(ni);
end

%y = cos(2*pi*x(:,1) + pi*x(:,2)');
%fexact = -2*sin(2*pi*x(:,1) + pi*x(:,2)');

y = cos(2*pi*x{1}) + sin(pi*x{2}');
fexact = -sin(2*pi*x{1}) + cos(pi*x{2}');

problem = "fbenney";
params = [1,1,1,1,1];

%problem = "fTest2";

method = "diff_ps_2d";

if method == "diff_fd"
    init_fd(x, 2);
end

tic
f = feval(problem,x,y,params,method);
tc = toc;

if method == "diff_fd"
    del_fd();
end



fe = f - fexact;

err = norm(fe,'inf');

fprintf('Problem: %s, Method: %s, Error: %g, Time: %gs\n',problem,method,err,tc);

[X,Y] = meshgrid(x{2},x{1});
surf(X,Y,real(f));
shading interp;

function f = fTest(x,y,method)
    deg = [1,2];
    dy = feval(method,x,y,deg);
    f = dy(:,2)/(4*pi^2) + dy(:,1)/(2*pi) + y;
end

function f = fTest2(x,y,method)
    deg = [1,0;0,1]';
    dy = feval(method,x,y,deg);
    f = dy(:,:,1)/2/pi + dy(:,:,2)/pi;
end
