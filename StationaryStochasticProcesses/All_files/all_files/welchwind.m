
function [welch,WEI]=welchwind(N,now);
% [welch,WEI]=welchwind(N,now);

 NN=ceil(2/(now+1)*N-1);
  step=(1/(now+1)*N);
  win=hanning(NN);
  win=win./sqrt(win'*win);
  welch=[win; zeros(N-NN,1)];
  if now>=2
    for i=2:now-1
      welch=[welch [zeros(fix(step*(i-1)),1) ;win; zeros(N-NN-fix(step*(i-1)),1)]];
    end
    welch=[welch [zeros(N-NN,1); win]];
  end  

   WEI=ones(now,1);



