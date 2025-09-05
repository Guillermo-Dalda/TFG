function val = artificial(x)
   val = (x(3) + 2.0) * x(1)*x(1) * x(2);
   g1 = 1 - ((x(2)*x(2)*x(2) * x(3)) / (71785 * x(1)*x(1)*x(1)*x(1)));
   g2 = (4 * x(2)*x(2) - x(1) * x(2)) / (12566 * (x(2) * x(1)*x(1)*x(1) - x(1)*x(1)*x(1)*x(1))) + (1 / (5108 * x(1)*x(1))) - 1;
   g3 = 1 - ((140.45 * x(1)) / (x(2)*x(2) * x(3)));
   g4 = (x(2) + x(1)) / 1.5 - 1;
   
   determinant = matrix_operations();

   val = determinant * val + g1 * (g1>0) + g2 * (g2>0) + g3 * (g3>0) + g4 * (g4>0);
end

function determinant = matrix_operations()
    N = 50;
    matrix = eye(N);
    aux = eye(N);
    determinant = 1;

    for i = 1:N
        for j = 1:N
            for k = 1:N
                matrix(i,k) = matrix(i,k) + aux(j,k) * j;
            end
        end
    end

    inverse = invert_matrix(matrix, N);
    
    aux = zeros(N);

    for i = 1:N
        for j = 1:N
            for k = 1:N
                aux(i,j) = aux(i,j) + matrix(i,k) * inverse(k,j);
            end
        end
    end
    
    for i = 1:N
        determinant = determinant * aux(i,i);
    end
end

function inverse = invert_matrix(matrix, N)
    inverse = eye(N);
    
    for i = 1:N
        pivot = matrix(i,i);
        for j = 1:N
            matrix(i,j) = matrix(i,j) / pivot;
            inverse(i,j) = inverse(i,j) / pivot;
        end
        for j = 1:N
            if i ~= j
                factor = matrix(j,i);
                for k = 1:N
                    matrix(j,k) = matrix(j,k) - matrix(i,k) * factor;
                    inverse(j,k) = inverse(j,k) - inverse(i,k) * factor;
                end 
            end
        end
    end
end