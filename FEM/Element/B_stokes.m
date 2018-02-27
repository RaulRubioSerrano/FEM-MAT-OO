classdef B_stokes < B
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods(Access = ?Element_stokes)
        function [obj,B] = computeB(obj,nunkn,nelem,nnode,cartd)
            
            B = zeros(2,nnode*nunkn,nelem);
            for i = 1:nnode
                j = nunkn*(i-1)+1;
                B(1,j,:)  = cartd(1,i,:);
                B(2,j+1,:)= cartd(1,i,:);
                B(3,j,:)  = cartd(2,i,:);
                B(4,j+1,:)= cartd(2,i,:);
            end
            obj.value = [obj.value {B}];
        end
    end
end

