% Return set of modes from samples sampled from these modes

function res = modes_rnd(samples)
	Tmodes = size(samples, 3);
	
	sample_rect = reshape(samples, 2, Tmodes)';

	% get some modes
	[uniq, ~, idx] = unique(sample_rect, "rows");
	cnts = accumarray(idx(:),1);
	umodes = [uniq cnts];
	staying_time = 8;
	modes = umodes(cnts > staying_time, :);
	lmodes = size(modes, 1);
	printf("Number of potential modes: %i\n", lmodes);
	if (lmodes == 1)
		error("A single mode, likely all steps have been rejected and no good starting point has been found, failure!\n");
	elif (lmodes < 5)
		lsmodes = size(umodes(cnts >= 2, :), 1);
		printf("There are %i locations with two or more rejection steps. If there are many of these locations the landscape is too flat and modes are difficult to find!\n");
	end

	% sort modes on number of cnts
	max_num_modes=10;
	ff = flipud(sortrows(modes,3));
	nmodes = min(lmodes, max_num_modes);
	topff = ff(1:nmodes, 1:2);

	%printf("The potential modes found (with colums [x y freq]):\n");
	res = topff;
end
