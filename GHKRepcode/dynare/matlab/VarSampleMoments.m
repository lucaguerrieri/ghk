function [YtY,XtY,YtX,XtX,Y,X] = VarSampleMoments(FirstObservation,LastObservation,qlag,trend);
global options_

X = [];
Y = [];
YtY = [];
YtX = [];
XtY = [];
XtX = [];


if exist(options_.datafile)
  eval(options_.datafile);
else
  eval(['load ' options_.datafile]);
end
data = [ ];
for i=1:size(options_.varobs,1)
    data = [data eval(deblank(options_.varobs(i,:)))];
end

if qlag > FirstObservation
    disp('VarSampleMoments :: not enough data to initialize! Try to increase FirstObservation.')
    return
end
NumberOfObservations = LastObservation-FirstObservation+1;
NumberOfVariables = size(data,2);
if trend == -1% No constant
    X = zeros(NumberOfObservations,NumberOfVariables*qlag);
elseif trend == 0;% Constant
    X = zeros(NumberOfObservations,NumberOfVariables*qlag+1);
    indx = NumberOfVariables*qlag+1;
elseif trend == 1;% Constant + Trend
    X = zeros(NumberOfObservations,NumberOfVariables*qlag+2);
    indx = NumberOfVariables*qlag+1:NumberOfVariables*qlag+2;
else
    disp('varols :: trend must be equal to -1,0 or 1!')
    return
end
% I build matrices Y and X  
Y = data(FirstObservation:LastObservation,:);
for t=1:NumberOfObservations
    line = t + FirstObservation-1;
    for lag = 1:qlag
        X(t,(lag-1)*NumberOfVariables+1:lag*NumberOfVariables) = data(line-lag,:);
    end
    if trend == 0
        X(t,indx) = 1;
    elseif trend == 1
        X(t,indx) = [ 1 , t ];
    end
end 

YtY = Y'*Y;
YtX = Y'*X;
XtY = X'*Y;
XtX = X'*X;