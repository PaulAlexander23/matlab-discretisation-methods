function init_fd(x, order)
    nx = length(x);
    dx = x(2)-x(1);
    global D
    if isempty(D)
        D = cell(1,4);
        if (order == 2)
            D{1} = spdiags(ones(nx,1)*[-1,1]/2,[-1,1],nx,nx)/dx;
            D{1} = spdiags(ones(nx,1)*[-1,1]/2,[-1,1] + nx,nx,nx)/dx+D{1};
            D{1} = spdiags(ones(nx,1)*[-1,1]/2,[-1,1] - nx,nx,nx)/dx+D{1};
            D{2} = spdiags(ones(nx,1)*[1,-2,1],[-1,0,1],nx,nx)/dx/dx;
            D{2} = spdiags(ones(nx,1)*[1,-2,1],[-1,0,1]+nx,nx,nx)/dx/dx+D{2};
            D{2} = spdiags(ones(nx,1)*[1,-2,1],[-1,0,1]-nx,nx,nx)/dx/dx+D{2};
            D{3} = spdiags(ones(nx,1)*[-1,2,0,-2,1]/2,[-2,-1,0,1,2],nx,nx)/dx/dx/dx;
            D{3} = spdiags(ones(nx,1)*[-1,2,0,-2,1]/2,[-2,-1,0,1,2]+nx,nx,nx)/dx/dx/dx+D{3};
            D{3} = spdiags(ones(nx,1)*[-1,2,0,-2,1]/2,[-2,-1,0,1,2]-nx,nx,nx)/dx/dx/dx+D{3};
            D{4} = spdiags(ones(nx,1)*[1,-4,6,-4,1],[-2,-1,0,1,2],nx,nx)/dx/dx/dx/dx;
            D{4} = spdiags(ones(nx,1)*[1,-4,6,-4,1],[-2,-1,0,1,2]+nx,nx,nx)/dx/dx/dx/dx+D{4};
            D{4} = spdiags(ones(nx,1)*[1,-4,6,-4,1],[-2,-1,0,1,2]-nx,nx,nx)/dx/dx/dx/dx+D{4};
        elseif (order == 4)
            D{1} =      spdiags(ones(nx,1)*[1,-8,0,8,-1]/12,[-2,-1,0,1,2],nx,nx)/dx;
            D{1} = D{1} + spdiags(ones(nx,1)*[1,-8,0,8,-1]/12,[-2,-1,0,1,2]+nx,nx,nx)/dx;
            D{1} = D{1} + spdiags(ones(nx,1)*[1,-8,0,8,-1]/12,[-2,-1,0,1,2]-nx,nx,nx)/dx;
            D{2} =      spdiags(ones(nx,1)*[-1,16,-30,16,-1]/12,[-2,-1,0,1,2],nx,nx)/dx/dx;
            D{2} = D{2} + spdiags(ones(nx,1)*[-1,16,-30,16,-1]/12,[-2,-1,0,1,2]+nx,nx,nx)/dx/dx;
            D{2} = D{2} + spdiags(ones(nx,1)*[-1,16,-30,16,-1]/12,[-2,-1,0,1,2]-nx,nx,nx)/dx/dx;
            D{3} =      spdiags(ones(nx,1)*[1,-8,13,0,-13,8,-1]/8,[-3,-2,-1,0,1,2,3],nx,nx)/dx/dx/dx;
            D{3} = D{3} + spdiags(ones(nx,1)*[1,-8,13,0,-13,8,-1]/8,[-3,-2,-1,0,1,2,3]+nx,nx,nx)/dx/dx/dx;
            D{3} = D{3} + spdiags(ones(nx,1)*[1,-8,13,0,-13,8,-1]/8,[-3,-2,-1,0,1,2,3]-nx,nx,nx)/dx/dx/dx;
            D{4} =      spdiags(ones(nx,1)*[-1,12,-39,56,-39,12,-1]/6,[-3,-2,-1,0,1,2,3],nx,nx)/dx/dx/dx/dx;
            D{4} = D{4} + spdiags(ones(nx,1)*[-1,12,-39,56,-39,12,-1]/6,[-3,-2,-1,0,1,2,3]+nx,nx,nx)/dx/dx/dx/dx;
            D{4} = D{4} + spdiags(ones(nx,1)*[-1,12,-39,56,-39,12,-1]/6,[-3,-2,-1,0,1,2,3]-nx,nx,nx)/dx/dx/dx/dx;
        end
    end
end
