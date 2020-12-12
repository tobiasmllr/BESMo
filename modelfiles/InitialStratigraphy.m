function [x_st,eta_st,pssi_st,pssi_sub] = InitialStratigraphy(LstrMat,nodes_N,MG,x,nLa,D90s,Msj,eta,delta_s,zbase,pssi0,psi)
format long e

%% Storing values of the stratigrahy in a matrix form. Active layer values are included. No structure is included within the substrate

x_st    = zeros(nodes_N,LstrMat);
eta_st  = zeros(nodes_N,LstrMat);
pssi_st = zeros(nodes_N,LstrMat,MG); 
pssi_sub = zeros(nodes_N,MG);

% Matrix of locations of the stratigrahy
for m = 1:LstrMat
    for n = 1:nodes_N
        x_st(n,m) = x(n);
    end
end

for n = 1:nodes_N
    La_str = nLa*D90s(n,1)/1000;
    for ii = 1:LstrMat        
        % below uppermost layer of the substrate---------------------------
        if ii < Msj(n,1)                                
            eta_st(n,ii) = delta_s*(ii - 1) + zbase;                                            
            for k = 1:MG                                
                pssi_st(n,ii,k) = pssi0(k);                                        
            end
        % uppermost layer of the substrate---------------------------------    
        elseif ii == Msj(n,1)  
            eta_st(n,ii) = eta(n,1)-La_str;            
            for k = 1:MG                    
                pssi_st(n,ii,k) = pssi0(k);                                        
            end                        
            for k = 1:MG                           
                pssi_sub(n,k) = pssi_st(n,ii,k);                
            end                        
        % Active layer-----------------------------------------------------                      
        elseif ii ==  Msj(n,1) + 1
            eta_st(n,ii) = eta(n,1);            
            for k = 1:MG                    
                pssi_st(n,ii,k) = psi(n,k);                                        
            end                       
        % data "above" the bed surface-------------------------------------
        else
            eta_st(n,ii) = NaN;            
            for k = 1:MG                    
                pssi_st(n,ii,k) = NaN;                                        
            end            
        end       
    end    
end
