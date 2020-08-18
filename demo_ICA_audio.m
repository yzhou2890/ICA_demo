function demo_ICA_audio()
    % demo of ICA
    % 
    % 
    
    % see Amari1996_demo() for the 3-source version
    
    % S     - uncorrupted source signal
    % X     - mixed signal, X = A*S, where A is the mixing matrix
    % Y     - demixed signal, Y = W*X, where W is the demixing matrix
    %                           = W*A*X

    % S     - d x n
    % X     - d x n
    % Y     - d x n
    
    
    % assumption:
    % Exp[y]    = 0
    % Exp[y^2]  = 1
    % num. of collected signals is no less than the num/dim of sources

global S X Y
global d Fs
    
    % video preparation
    hvd = VideoWriter(['X_',num2str(round(rand*1e5)),'.avi']);
    hvd.FrameRate = 3;
    hvd.Quality = 95;
%    open(hvd);
    
 
    N = 1e4;
    load chirp;  s1 = y(1:N); s1 = s1/std(s1);
    load gong;   s2 = y(1:N); s2 = s2/std(s2);
    S = [s1, s2, rand(N,1)]';
    d = size(S,1);
    n = size(S,2);
    A = rand(d);
    Fs = 1e4;
    X = A*S;    
    
                  

    str = {{'signal-1','signal-2','signal-3'},{'1st mix','2nd mix', '3rd mix'},...
           {'1st demix','2nd demix','3rd demix'}};
    
    
    
    
    % prepare to show iterations
    figure(2), 
    set(gcf,'color','w','Position',[100 50 1300, 195*d]), 
    set(gcf,'WindowButtonDownFcn',{@MyWindowButtonDownFcn});
    set(gcf,'menubar','None');
    
    lft = linspace(0.03,1,4); lft(end) = []; 
    btm = linspace(0.04,1,d+1); btm(end) = [];  btm = fliplr(btm); 
    for i=1:3
        for k =1:d                                    
            s1 = repmat([ones(10),zeros(10); zeros(10), ones(10)],5,11);
            subplot(d,3,3*k-2, 'Position',[lft(1),btm(k), 1/3*0.9, 1/d*0.82 ]),
            imshow( s1,[-5,2]);
            %text( 40,50,'?','background','w','FontSize',80)
            title(str{1}{k})

            s1 = X(k,:);                
            subplot(d,3,3*k-1, 'Position',[lft(2),btm(k), 1/3*0.9, 1/d*0.82 ]),
            plot(s1,'r'), ylim([min(s1)-0.1*abs(min(s1)), max(s1)+0.1*abs(max(s1))]*1);
            legend(str{2}{k})
            
            s1 = repmat([ones(10),zeros(10); zeros(10), ones(10)],5,11);
            subplot(d,3,3*k,   'Position',[lft(3),btm(k), 1/3*0.9, 1/d*0.82 ]),
            imshow( s1,[-5,2]);
            %text( 0, 0,'','background','w')            
            title(str{3}{k})
        end
%        writeVideo(hvd, getframe(gcf));
    end
    
    % there are two options:
    % 1:    SIGNAL SEPARATION VIA [AMARI 1996], key equation: (19)
    % 2:    [Haykin 2001, p510-p525], key-equation: (10.104) in p521
    %       implement batch mode of (10.104) [see P544 of Haykin01]       
    
    eta = 0.1;                  % learning rate
    W0  = rand(d)*0.02 - 0.01;  % initial demxing matrix

    for i=1:200
        if 1
            Y    = W0*X ;           % d x n                        
            PhiY = zeros(d, n);     % d x n

            for j = 1:d
                for k =1:n                                                
                    %PhiY(j,k) =   Haykin01_Eqn_10p93( Y(j,k) ) ;   % Haykin2001,p521,eqn 10.93
                    PhiY(j,k) =   Amari96_Eqn_14( Y(j,k) ) ;        % Amari96

%                    W0 = W0 + eta * (eye(d) -  PhiY(j,k) * Y(j,k)) * W0;  % stohcastic version
                end
%                W0 = W0 + eta * (eye(d) -  PhiY(j,:) * Y(j,:)') * W0;  % stohcastic version
            end
        end
        W0 = W0 + eta * (eye(d) - 1/n * PhiY * Y') * W0;     % batch version
        

        if 0 % fast-ica
            addpath('FastICA_25')
            [Out1, A_, W_] = fastica(X);

            Y  = Out1;
            W0 = W_;
            disp( A_ ),
            disp( W_ ),        
        end            


        % show mid-term results to generate a video for demo
        if rem(i,5)==0        
            figure(2)
            set(gcf,'name',[num2str(i),'th iteration:'])               

            for k =1:d                
                s1 = Y(k,:);                
                subplot(d,3,3*k,   'Position',[lft(3),btm(k), 1/3*0.9, 1/d*0.82 ]), 
                plot(s1,'b'), ylim([min(s1)-0.1*abs(min(s1)), max(s1)+0.1*abs(max(s1))]*1);
                legend([num2str(i),'th iteration - ', str{3}{k}],'background','w')
            end

%            writeVideo(hvd, getframe(gcf));
        end
    end        
    
    
    for i=1:3
        for k =1:d
            s1 = S(k,:);                                       
            subplot(d,3,3*k-2, 'Position',[lft(1),btm(k), 1/3*0.9, 1/d*0.82 ]),
            plot(s1,'g'), ylim([min(s1)-0.1*abs(min(s1)), max(s1)+0.1*abs(max(s1))]*1);       
            legend(str{1}{k})
        end
 %       writeVideo(hvd, getframe(gcf));
    end    
 %   close(hvd);

    
    % corre. in-betwtn souces and ICA's outputs          
    mCrr = zeros(size(S,1), size(Y,1));
    for i=1:size(S,1)
        for j=1:size(Y,1)
            mCrr(i,j) = sum(S(i,:).*Y(j,:))/(norm(S(i,:))*norm(Y(j,:)));
        end
    end
    disp('Correlation: (S_i, Y_j)')
    disp(mCrr);
    
    
    % error indicator (Eqn. in Sec. 5 of [Amari 96]) for the system
    P  = abs( W0 * A ); 
    P1 = P/(diag(max(P,[],1)));
    P2 = (diag(max(P,[],2)))\P;
    E = sum(P1(:)) + sum(P2(:)) -  size(P,1) - size(P,2)
 
        
    
    % DISPLAY RESULTS
    disp('A = ')
    disp(A)
    disp('W0 = ')
    disp(W0)
    disp('W0*A = ')
    disp(W0*A)
    
    % 
    for k = 1:d
 %       soundsc(Y(k,:),Fs)
    end

end


function fy = Amari96_Eqn_14(y)
%  Eqn (14) in [Amari96]
    vcoeff = [ 3/4,  15/4,  -14/3,  -29/4,  29/4 ];
%    vcoeff = [ 3/4, 25/4,  -14/3,  -47/4,  29/4 ];  % exact coef. in Eqn (14) of [Amari96]
    vy     = [ y^11,  y^9,    y^7,    y^5,   y^3 ];
    fy     = sum( vcoeff.*vy );
end


function fy = Haykin01_Eqn_10p93(y)
% Eqn (10.93) in [Haykin 2001] at page 519
    vcoeff = [ 1/2, 2/3, 15/2, 2/15, -112/3,  128, -512/3];
    vy     = [ y^5, y^7,  y^9, y^11,   y^13, y^15,   y^17];
    fy     = sum( vcoeff.*vy );
end


function  MyWindowButtonDownFcn(~,~)
global S X Y    
global d Fs

    box = get(gcf,'Position');
    v = get(gcf,'CurrentPoint');
    
    i = round( v(1) * (3 - 1) / box(3)) + 1;
    j = d - round( v(2) * (d - 1) / box(4)) ;
    
    switch i
        case 1
            soundsc(S(j,:),Fs)
        case 2
            soundsc(X(j,:),Fs);
        case 3
            soundsc(Y(j,:),Fs);
    end
end
