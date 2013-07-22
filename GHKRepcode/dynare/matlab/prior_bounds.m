function bounds = prior_bounds(bayestopt)
  global options_
  
  pshape = bayestopt.pshape;
  pmean = bayestopt.pmean;
  p1 = bayestopt.p1;
  p2 = bayestopt.p2;
  p3 = bayestopt.p3;
  p4 = bayestopt.p4;
  
  n = length(pmean);
  bounds = zeros(n,2);
  
  for i=1:n
    switch pshape(i)
     case 1
      mu = (pmean(i)-p3(i))/(p4(i)-p3(i));
      stdd = p2(i)/(p4(i)-p3(i));
      A = (1-mu)*mu^2/stdd^2 - mu;
      B = A*(1/mu - 1);
      bounds(i,1) = qbeta(options_.prior_trunc,A,B)*(p4(i)-p3(i))+p3(i);
      bounds(i,2) = qbeta(1-options_.prior_trunc,A,B)*(p4(i)-p3(i))+p3(i);
     case 2
      b = p2(i)^2/(pmean(i)-p3(i));
      a = (pmean(i)-p3(i))/b;
      bounds(i,1) = mj_qgamma(options_.prior_trunc,a)*b+p3(i);
      bounds(i,2) = mj_qgamma(1-options_.prior_trunc,a)*b+p3(i);
     case 3
      bounds(i,1) = qnorm(options_.prior_trunc,pmean(i),p2(i));
      bounds(i,2) = qnorm(1-options_.prior_trunc,pmean(i),p2(i));
     case 4
      nu = p2(i);
      mu = pmean(i);
      beta = ( gamma( (nu-1)/2 ) / mu / gamma( nu/2 ) )^2;
      a=2/beta;
      bounds(i,1) = 1/sqrt(mj_qgamma(1-options_.prior_trunc,p2(i)/2)*beta);
      bounds(i,2) = 1/sqrt(mj_qgamma(options_.prior_trunc,p2(i)/2)*beta);
     case 5
      bounds(i,1) = p1(i);
      bounds(i,2) = p2(i);
     otherwise
      bounds(i,1) = -Inf;
      bounds(i,2) = Inf;
    end
  end
  