function n = unifpdfN(X, A, B)
	N=length(X);

	for i = 1:N
		a=A(i);
		b=B(i);
		x=X(i);
		if (a > b)
			[b,a] = deal(a,b);
		end
		n(i) = prod(unifpdf(x, a, b));
	end
	n=n';
end
