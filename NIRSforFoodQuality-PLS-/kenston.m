							
%#  function [object,xm,ym,xt,yt]=kenston(x,no_p,men,pl,y)				
%#									
%#  AIM: 	Kennard-Stone design (e.g. for subset selection).	
%#									
%#  PRINCIPLE:  Based on Computer Aided Design of Experiments. 		
%#		The first point can be the mean or the furthest from 	
%#		the mean.						
%# 		REF :   R. W. Kennard and L. A. Stone			
%# 			Technometrics Vol. 11, No. 1, 1969 		
%# 									
%#  INPUT:	x :  (n x m) absorbent matrix with n spectra		
%#			 and m variables				
%#		no_p : number of objects to be selected			
%#		men: position of the first point selected		
%#		     (1 closest to mean; 0 furthest from mean)		
%#		pl: 1 plot; 0 no plot					
%#		y: matrix of responses							
%#  OUTPUT:	object : (1 x no_p) the vector of designed objects	
%#																	
%#  AUTHOR: 	Wen Wu 							
%#	    	Copyright(c) 1997 for ChemoAC				
%#          	FABI, Vrije Universiteit Brussel            		
%#          	Laarbeeklaan 103 1090 Jette				
%#    	    								
%# VERSION: 1.1 (28/02/1998)						
%#									
%#  TEST:   	Roy De Maesschalck					
%#									

function [object,xm,ym,xt,yt]=kenston(x,no_p,men,pl,y);

[n,m]=size(x);	
t=x;
% Kennard and Stone method to select objects
  % starting (1st) point is the closest (or furthest) from the centroid
	meant=mean(t);
	t1=t-ones(n,1)*meant;
    for i=1:n
        a(i)=t1(i,:)*t1(i,:)'; 
	end %i
	if men==1,
		[b,c]=min(a);
	else
		[b,c]=max(a);
	end
	object(1)=c;
	clear a b c t1
   % 2nd points is the furthest point from 1st point
	t1=t-ones(n,1)*t(object(1),:);
	for i=1:n
	    a(i)=t1(i,:)*t1(i,:)';
	end %i
	[b,c]=max(a);
	object(2)=c;
	clear a b c t1	
  % k+1 point
	for pi=3:no_p
	    list=1:n;
	    tt=t;
	    k=length(object);
	    list(object)=[];
	    tt(object,:)=[];
	    nl=length(list);
	    for j=1:nl
		for i=1:k
	 	    t1=tt(j,:)-t(object(i),:);
		    a(i)=t1*t1';
		end % i
	    	[b,c]=min(a);
		dmin(j)=b;
	    end %j
	    [b,c]=max(dmin');
	    object(pi)=list(c);
	    clear dmin a b c list
	end %pi
object=object(1:no_p);
% plot
if pl==1
	[a,b]=size(t);
   if b>1
	plot(t(:,1),t(:,2),'.')
	for i=1:n
	    text(t(i,1),t(i,2),int2str(i))
	end
	hold on
	plot(t(object,1),t(object,2),'r*')
	hold off
	xlabel('Variable 1')
	ylabel('Variable 2')
   end 
   if b==1
	plot(t(:,1),t(:,1),'.')
	for i=1:n
	    text(t(i,1),t(i,1),int2str(i))
	end
	hold on
	plot(t(object,1),t(object,1),'r*')
	hold off
	xlabel('Variable 1')
	ylabel('Variable 1')
   end 
end                    
xm=x(object,:);
ym=y(object,:);
ind=[1:n]';
ind(object)=[];
xt=x(ind,:);
yt=y(ind,:);