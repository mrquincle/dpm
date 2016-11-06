function disp_reverse(nd_object, parent_object = [])
	nd_object
	ndim = ndims(nd_object)
	if (ndim > 2)
		nsize = size(nd_object,1)
		for i = 1:nsize
%			new_object = nd_object(:,:,i)
			new_object = nd_object(i,:,:,:,:)
			new_object = nd_object(i)
			error("test")
			disp_reverse(new_object, [parent_object i]);
		end
	else
		parent_object = nd_object
	end
end
