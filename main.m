
n = 2^3;
x = linspace(1/n,1,n);

y = cos(2*pi*x)';

problem = @(x,y,dy) testEquation(y,dy(:,1),dy(:,2));
problemDeg = [1,2];

method = "diff_fd";

if method == "diff_fd"
    del_fd()
    init_fd(x, 2);
end

tic
f = problem(x,y,feval(method,x,y,problemDeg));
tc = toc;

if method == "diff_fd"
    del_fd();
end

fexact = -sin(2*pi*x)';

fe = f - fexact;

err = norm(fe,'inf');

fprintf('Problem: %s, Method: %s, Error: %g, Time: %gs\n',func2str(problem),method,err,tc);

plot(x,f,x,fexact);


function f = testEquation(y,dy,d2y)
    f = d2y/(4*pi^2) + dy/(2*pi) + y;
end
