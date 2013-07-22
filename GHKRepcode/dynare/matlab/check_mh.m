function check_mh(fname)
  eval(['load ' fname]);
  nb = size(x3,3);
  np = size(x2,2);
  nr = size(x2,1);
  
  j1 = ceil(0.5*nr);
  x = [j1:100:nr];
  z = [];
  for i=1:np
    y1 = zeros(size(x),nb);
    for k=1:nb
      for j=1:length(x)
	y1(j,k) = mean(x2(1:x(j),i,k));
      end
    end
    my_subplot(i,np,4,5,'MH convergence');
    plot([y1])
    xmin = min(min(x2(:,i,:)));
    xmax = max(max(x2(:,i,:)));
    z = [z; [i sum(sum(x3(:,i,:) < xmin))/nr sum(sum(x3(:,i) > xmax))/nr]];
  end
  disp(z)
     