function D = init_fd2(x, degree, accuracy)
    dim = size(degree,1);
    if dim > 1
        nx = cellfun(@(x) length(x),x);
        dx = cellfun(@(x) x(2)-x(1),x);
    else
        nx = length(x);
        dx = x(2)-x(1);
    end
    
    D = cell(1,size(degree,2));
    for j = 1:size(degree,2)
        D{j} = 1;
        for k = 1:dim
            [s,c] = centred_scheme(degree(k,j), accuracy);
            
            Dmat = construct_matrix(nx(k),dx(k),s,c,degree(k,j));
            
            D{j} = kron(Dmat,D{j});
        end
    end
    
    function [s,c] = centred_scheme(degree, accuracy)
        sL = ceil(degree/2) + ceil(accuracy/2)-1;
        s = -sL:sL;
        ex = (0:2*sL)';
        S = repmat(s,sL*2+1,1).^ex;
        del = factorial(degree) * (ex == degree);
        c = (S\del)';
    end
    
    function mat = construct_matrix(nx,dx,s,c,d)
        mat = spdiags(ones(nx,1)*c,s,nx,nx)*dx.^(-d);
        mat = mat + spdiags(ones(nx,1)*c,s+nx,nx,nx)*dx.^(-d);
        mat = mat + spdiags(ones(nx,1)*c,s-nx,nx,nx)*dx.^(-d);
    end
end