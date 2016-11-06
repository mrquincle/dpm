% Construct pdf from differences between the modes

function res = mode_proposal_pdf(modes)

	% get diffs between modes as jump proposal distribution
	[tm tn] = size(modes);
	pdf_diff = reshape(bsxfun(@minus,reshape(modes,[1 tm tn]),reshape(modes,[tm 1 tn])),[tm*tm tn]);

	% we need can decide to interpolate between modes as well
	%pdf_diff = [pdf_diff ; pdf_diff / 2];

	% TODO: This step has an arbitrary parameter value 5, let's see if we actually need this?
	% get rid of jumps that are too large
	smalljumps = abs(pdf_diff(:,1) .* pdf_diff(:,2)) < 5;
	pdf_diff = pdf_diff(find(smalljumps),:);

	% get rid of jumps that are no jumps (difference is 0)
	smalljumps = abs(pdf_diff(:,1) .* pdf_diff(:,2)) != 0;
	pdf_diff = pdf_diff(find(smalljumps),:);

	% TODO: This step has an arbitrary parameter value 5, let's see if we actually need this?
	if(size(pdf_diff,1) < 5)
		error("No success in finding modes! Run this piece of code again, perhaps let the chains walk for a bit longer.\n");
	end

	res = pdf_diff;
end
