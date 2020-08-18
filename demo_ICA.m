function demo_ICA()
    % 1st demo for ICA
    % 
    % 
    
    % S     - uncorrupted source signal
    % X     - mixed signal, X = A*S, where A is the mixing matrix
    % Y     - demixed signal, Y = W*X, where W is the demixing matrix
    %                           = W*A*X

    % S     - 3 x n
    % X     - 3 x n
    % Y     - 3 x n
    
    
    % assumption:
    % Exp[y]    = 0
    % Exp[y^2]  = 1
    % num. of collected signals is no less than the num/dim of sources
    

    hvd = VideoWriter(['X_',num2str(round(rand*1e5)),'.avi']);
    hvd.FrameRate = 10;
    hvd.Quality = 75;
%    open(hvd);
    
    
    % SIGNAL GENERATION
    t = (0:0.0001:0.5);     % time 
    n = length(t);          % number of samp. points
    
    % source signals - {s1, s2, s3}
%    S = [rand(1,n)*2 - 1;  0.1*sin(400*t).*cos(30*t); 0.01*sign( sin( 500*t + 9*cos(40*t) )) ];
    S = [randn(1,n)*2 - 1;  1e-1*sin(400*t).*cos(30*t); 1e-2*sign( sin( 500*t + 9*cos(40*t) )) ]*1;    
%    S = [rand(1,n)*2 - 10;  0.1*sin(400*t).*cos(30*t)+0; 0.01*sign( sin( 500*t + 9*cos(40*t) )); exp(-40*t).*(1+sin(100*t)) ];


    

    
    % SIGNAL MIXING matrix - A
    d = size(S,1);
    A = rand(d,d)*0.02     % mixing matrix
%    A = [ 0.56 0.79 -0.37; -0.75 0.65 0.86; 0.17 0.32 -0.48 ]; % [Haykin01,p523]
%     A = [   0.0016    0.0146    0.0004
%             0.0027    0.0142    0.0181
%             0.0171    0.0107    0.0079 ]+1;   % difficult mixing matrix
%    A = eye(d) + 0.01*rand(d);


    X = A*S;                % mixed signals, i.e., raw signals from sensors
%    tmp = X - repmat(mean(X,2), 1, n);
%    X = inv( diag(sqrt( sum(tmp.*tmp, 2)/(n-1) )) )* tmp;


%     figure(1), set(gcf,'color','w'),
%     subplot(3,2,1), plot(t,s1), legend('s_1'), ylim([min(s1)-0.1*abs(min(s1)), max(s1)+0.1*abs(max(s1))]*1);
%     subplot(3,2,3), plot(t,s2), legend('s_2'), ylim([min(s2)-0.1*abs(min(s2)), max(s2)+0.1*abs(max(s2))]*1);
%     subplot(3,2,5), plot(t,s3), legend('s_3'), ylim([min(s3)-0.1*abs(min(s3)), max(s3)+0.1*abs(max(s3))]*1);
% 
%     subplot(3,2,2), plot(t,X(1,:)), legend('x_1'), ylim([min(X(1,:))-0.1*abs(min(X(1,:))), max(X(1,:))+0.1*abs(max(X(1,:)))]*1);
%     subplot(3,2,4), plot(t,X(2,:)), legend('x_2'), ylim([min(X(2,:))-0.1*abs(min(X(2,:))), max(X(2,:))+0.1*abs(max(X(2,:)))]*1);
%     subplot(3,2,6), plot(t,X(3,:)), legend('x_3'), ylim([min(X(3,:))-0.1*abs(min(X(3,:))), max(X(3,:))+0.1*abs(max(X(3,:)))]*1);

    
    % prepare online illustration
    figure(2), set(gcf,'color','w','position',[100, 50, 1000, 600]),
    for k=1:d 
        s1 = X(k,:);
        subplot(d,2,2*k-1), plot(t,s1,'r'), legend(['x_',num2str(k)]), 
        ylim([min(s1)-0.1*abs(min(s1)), max(s1)+0.1*abs(max(s1))]*1);
    end
    
   
    % there are two options:
    % 1:    SIGNAL SEPARATION VIA [AMARI 1996], key equation: (19)
    % 2:    [Haykin 2001, p510-p525], key-equation: (10.104) in p521
    %       implement batch mode of (10.104) [see P544 of Haykin01]       
    
    if 1 
        eta = 0.1;                  % learning rate
        
%        W0  = eye(d);                % initial demixing matrix
        W0  = rand(d)*0.2 - 0.1;  % initial demxing matrix

        for i=1:300
            Y    = W0*X ;           % 3 x n                        
            PhiY = zeros(d, n);     % 3 x n
            
%              tmp = Y - repmat( mean(Y,2), 1, n);
%              Y = 2*diag(1./(sqrt(diag(cov(tmp')))))*tmp;
            
            for j = 1:d
                for k =1:n                                                
                    %PhiY(j,k) =   Haykin01_Eqn_10p93( Y(j,k) ) ;   % Haykin2001,p521,eqn 10.93
                    PhiY(j,k) =   Amari96_Eqn_14_with_4th_order_of_Eqn3( Y(j,k) ) ;    
                    %PhiY(j,k) =   Amari96_Eqn_14_with_5th_order_of_Eqn3(Y(j,k));
                    %PhiY(j,k) =   Bartlett02_Eqn_3(Y(j,k));
                    
%                    W0 = W0 + eta * (eye(d) -  PhiY(j,k) * Y(j,k)/n) * W0; %stochastic version, see below for batch version
                end                
            end           
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             

%             W1 = W0 + eta   *(eye(d) - 1/n * PhiY * Y') * W0;    % batch version here, see above for stochastic version 
%             disp(norm((eye(d) - 1/n * PhiY * Y') * W0)),
%             W0 = W1;
            
            W0 = W0 + eta   *(eye(d) - 1 * PhiY * Y' / norm(PhiY * Y')) * W0;    % batch version here, see above for stochastic version 
            W0 = W0 / norm( W0 );            
            
            % SHOW RESULTS ONLINE
            Y = W0*X;
            figure(2), set(gcf,'color','w'),
            k=1; s1 = Y(k,:); subplot(3,2,2), plot(t,s1), legend('y_1'), 
            title([num2str(i),'th iteration:'],'background','w'); % ylim( [min(-2,min(s1)), max(2,max(s1))]);
            k=2; s1 = Y(k,:); subplot(3,2,4), plot(t,s1), legend('y_2'), 
            title([num2str(i),'th iteration:'],'background','w'); % ylim( [min(-2,min(s1)), max(2,max(s1))]);
            k=3; s1 = Y(k,:); subplot(3,2,6), plot(t,s1), legend('y_3'), 
            title([num2str(i),'th iteration:'],'background','w'); % ylim( [min(-2,min(s1)), max(2,max(s1))]);            
            
%            writeVideo(hvd, getframe(gcf));
        end                
    end
    
    
    if 0 
        % fast-ica
        addpath('../ICA_Finland_FastICA_25/')
        [Out1, A_, W_] = fastica(X);
        
        Y  = Out1;
        W0 = W_;
        disp( A_ ),
        disp( W_ ),        
    end

    
    % error indicator (Eqn. in Sec. 5 of [Amari 96]) for the system
    P  = abs( W0 * A ); 
    P1 = P/(diag(max(P,[],1)));
    P2 = (diag(max(P,[],2)))\P;
    E = sum(P1(:)) + sum(P2(:)) -  size(P,1) - size(P,2);
        
        
    % SHOW RESULTS 
    figure(2), set(gcf,'color','w'),
    for k =1:d
        s1 = S(k,:);
        subplot(d,2,2*k-1), plot(t,s1,'g'), legend(['s_',num2str(k)]), 
        ylim([min(s1)-0.1*abs(min(s1)), max(s1)+0.1*abs(max(s1))]*1);
    end     
    
    for k = 1:10
 %       writeVideo(hvd, getframe(gcf));
    end
 %   close(hvd);
    
    % DISPLAY RESULTS
    disp('A = ')
    disp(A)
    disp('W0 = ')
    disp(W0)
    disp('W0*A = ')
    disp(W0*A)

end


function fy = Amari96_Eqn_14_with_4th_order_of_Eqn3(y)
%  Eqn (14) in [Amari96]
    vcoeff = [ 3/4,  15/4,  -14/3,  -29/4,  29/4 ];       
    pn     =  11:-2:3 ;
    fy     = sum( vcoeff.*(y.^pn) );
end


