function [ DtJtil ] = cal_DtJtil_min_inv_stiff(smoothDtC,matCh,DtC,dim,ngaus,alpha,hnorm,Msmooth,coordinatesn,coordinatesa,element,problembsc)
% 

npnod=dim.npnod; ndime=dim.ndime; nelem=dim.nelem; nnode=dim.nnode; 

strain0 = alpha;
DtJtil = [];

if (smoothDtC==1)
    [DtC] = smooth_DtC(DtC,3,npnod,ngaus,nelem,ndime,nnode,Msmooth,element,coordinatesn,coordinatesa,problembsc);

    DtC1 = zeros(npnod,1);
    for a=1:3
        for b=1:3
            DtC1(:) = DtC(a,b,:);
            DtJtil = DtJtil - strain0(a)*DtC1(:)*strain0(b)/((strain0*matCh*strain0')^2)/abs(hnorm)*(strain0*strain0');
        end
    end
    
elseif (smoothDtC==0)
    DtC1 = zeros(ngaus,nelem);
    for a=1:3
        for b=1:3
            DtC1(:,:) = DtC(a,b,:,:);
            DtJtil = DtJtil - strain0(a)*DtC1*strain0(b)/((strain0*matCh*strain0')^2)/abs(hnorm)*(strain0*strain0');
        end
    end
    
end

end

