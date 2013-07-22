a = [0.95 0.025; 0.025 0.95];

ngrid = 512;
freqs = 0 : ((2*pi)/ngrid) : (2*pi*(1 - .5/ngrid)); 

im = sqrt(-1);
  mathp_col = [];

    for ig = 1:ngrid,
      tpos  = exp( im*freqs(ig));
      tneg  =  exp(-im*freqs(ig));
      f_omega  =(1/(2*pi))*[inv(eye(2)-a*tneg);eye(2)]* ...
		Sigma_e_*[inv(eye(2)-a'*tpos) eye(2)];
      f_hp = f_omega;
      mathp_col = [mathp_col ; (f_hp(:))'];    % store as matrix row
                                               % for ifft
    end;

    % covariance of filtered series
    imathp_col = real(ifft(mathp_col))*2*pi;
    reshape(imathp_col(1,:),4,4)
