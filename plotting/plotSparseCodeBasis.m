function plotSparseCodeBasis(global_info, basis_vectors)
n = numel(basis_vectors);
for i = 1:n
   plotSparseCodeAtom(global_info, 0, basis_vectors(i).direction_x, basis_vectors(i).direction_y); 
end


end