% -- Function File: T = meshtable (X, Y, Z)
%     Reconstruct data in table form. 
%
%     The input for surf or mesh is the same as that of meshtable. It will be
%     casting X and Y into [X Y] pairs in all possible combinations, expecting
%     Z to be ordered accordingly.
%
function T = meshtable(X, Y, Z)
	[XX YY] = meshgrid(X, Y);
	T = [XX(:) YY(:) Z(:)];
end
