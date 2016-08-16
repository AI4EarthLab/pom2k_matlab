function [h] = slpmax(h,im,jm,fsm,slmax)
% **********************************************************************
% *                                                                    *
% * FUNCTION    :  Limits the maximum of:                              *
% *                                                                    *
% *                  <difference of depths>/<sum of depths>            *
% *                                                                    *
% *                for two adjacent cells. The maximum possible value  *
% *                is unity.                                           *
% *                                                                    *
% **********************************************************************

for loop=1:10
%     Sweep right:
	for j=2:jm-1
		for i=2:im-1
        	if(fsm(i,j) ~= 0.e0 && fsm(i+1,j) ~=0.e0) 
        		if(abs(h(i+1,j)-h(i,j))/(h(i,j)+h(i+1,j)) >= slmax)
            		mean=(h(i+1,j)+h(i,j))/2.e0;
                	%del=sign(slmax,h(i+1,j)-h(i,j));
                    if(h(i+1,j)-h(i,j)>=0.e0)
                        del = abs(slmax);
                    else del = -abs(slmax);
                    end
                	h(i+1,j)=mean*(1.e0+del);
                	h(i,j)=mean*(1.e0-del);
            	end
        	end
     	end 

%    Sweep left:
		for i=im-1:-1:2
    		if(fsm(i,j) ~= 0.0 && fsm(i+1,j) ~= 0.0) 
    			if(abs(h(i+1,j)-h(i,j))/(h(i,j)+h(i+1,j)) >= slmax)
            		mean=(h(i+1,j)+h(i,j))/2.e0;
            		%del=sign(slmax,h(i+1,j)-h(i,j));
                    if(h(i+1,j)-h(i,j)>=0.e0)
                        del = abs(slmax);
                    else del = -abs(slmax);
                    end
					h(i+1,j)=mean*(1.e0+del);
					h(i,j)=mean*(1.e0-del);
				end
    		end
    	end
	end
 
%   Sweep up:
	for i=2:im-1
    	for j=2:jm-1
        	if(fsm(i,j) ~= 0.0 && fsm(i,j+1) ~= 0.0) 
        		if(abs(h(i,j+1)-h(i,j))/(h(i,j)+h(i,j+1)) >= slmax) 
            		mean=(h(i,j+1)+h(i,j))/2.e0;
                    if(h(i,j+1)-h(i,j)>=0.e0) 
                        del=abs(slmax);
                    else del=-abs(slmax);
                    end 
                	%del=sign(slmax,h(i,j+1)-h(i,j));
                	h(i,j+1)=mean*(1.0+del);
                	h(i,j)=mean*(1.0-del);
            	end
        	end
		end
 
%   Sweep down:
		for j=jm-1:-1:2
			if(fsm(i,j) ~= 0.0 && fsm(i,j+1) ~= 0.0)
        		if(abs(h(i,j+1)-h(i,j))/(h(i,j)+h(i,j+1)) >= slmax)
            		mean=(h(i,j+1)+h(i,j))/2.e0;
                    if(h(i,j+1)-h(i,j)>=0.e0) 
                        del=abs(slmax);
                    else del=-abs(slmax);
                    end 
                	%del=sign(slmax,h(i,j+1)-h(i,j));
                	h(i,j+1)=mean*(1.e0+del);
                	h(i,j)=mean*(1.e0-del);
            	end
        	end
    	end
	end
end
