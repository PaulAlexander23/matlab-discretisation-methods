function init_fd_2d(x, order)
    n = cellfun(@(x) length(x),x);
    d = cellfun(@(x) x(2)-x(1),x);
    global D
    
    D = cell(1,4);
    if (order == 2)
        coeff = cell(5,1);
        coeff{1} = [-1,0,1]/2;
        coeff{2} = [-1,0,1]/2;
        coeff{3} = [1,-2,1];
        coeff{4} = [-1,2,0,-2,1]/2;
        coeff{5} = [1,-4,6,-4,1];
        
        for j = 1:5
            for k = 1:5
                Dx = diffmat(n(1),d(1),coeff{j},j-1);
                Dy = diffmat(n(2),d(2),coeff{k},k-1);
                
                D{j,k} = kron(Dx,Dy);
                
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
    
    function dmat = diffmat(nx,dx,coeff,order)
        diags = 1:length(coeff)-(length(coeff)+1)/2;
        dmat = spdiags(ones(nx,1)*coeff,diags,nx,nx)/(dx^order);
        dmat = dmat + spdiags(ones(nx,1)*coeff,diags + nx,nx,nx)/(dx^order);
        dmat = dmat + spdiags(ones(nx,1)*coeff,diags - nx,nx,nx)/(dx^order);
    end
end
