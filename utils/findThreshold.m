function [th] = findThreshold(shape, areas,partArea, w)

%peforms binary search
th = 0.5;
thmax = 1;
thmin = 0.0;
done = 0;
while(done == 0)
	tw = w;
	tw(tw <= th) = 0;
	tw(tw ~= 0) = 1;
	if(sum(tw .* areas) <=	partArea*0.95)
		%threshold is to high
		newTh = (th + thmin)/2;
		thmax = th;
		if(abs(newTh - th) < 0.001)
			done = 1;
		end
		th = newTh;
		continue
	end
	if(sum(tw .* areas) >	partArea)
		%threshold is to low
		newTh = (th + thmax)/2;
		thmin = th;
		if(abs(newTh - th) < 0.001)
			done = 1;
		end
		th = newTh;
		continue
	end
	%threshold is ok
	done = 1;
end



end
