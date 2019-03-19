function init_fd_2d(x, order)
    n = cellfun(@(x) length(x),x);
    dx = cellfun(@(x) x(2)-x(1),x);
    global D
    
    D = cell(1,4);
    if (order == 2)
        coeff = cell(4,1);
        coeff{1} = [-1,0,1]/2;
        coeff{2} = [1,-2,1];
        coeff{3} = [-1,2,0,-2,1]/2;
        coeff{4} = [1,-4,6,-4,1];
        
        for j = 1:4
            for iaxis = 1:2
                diags = 1:length(coeff{j})-(length(coeff{j})+1)/2;
                nx = n(iaxis);
                dx = d(iaxis);
                temp = spdiags(ones(nx,1)*coeff{j},diags,nx,nx)/dx;
                temp = temp + spdiags(ones(nx,1)*coeff{j},diags + nx,nx,nx)/dx;
                temp = temp + spdiags(ones(nx,1)*coeff{j},diags - nx,nx,nx)/dx;
                if iaxis == 1
                    D{j,iaxis} = kron(eye(n(2)),temp);
                else
                    D{j,iaxis} = kron(temp,eye(n(1)));
                end
            end
        end
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
