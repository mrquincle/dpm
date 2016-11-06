function show_schedule(data, k)

	palette = jet(k);

	y = data;
	x = [0:23];
	xi = [0:0.01:23];
	yi = interp1(x, y, xi, "cubic");
	area(xi, yi, 'facecolor', palette(k,:));

end


