function addcohesive_2d(NodeDict,ElemDict,ElsetPartsArray,ElsetMatrix)
#=
The function addcohesive_2d return the new cohesive elements in the dicts and the mapping dictionary between old and new nodes.
=#
    
    return 
end

function findsharednodes(ElemDict,Elset1,Elset2)
#=
The function return the nodes shared by the inclusions Elset1 and matrix Elset2.
=#  
    Nset1 = obtainnodes(ElemDict,Elset1);
    Nset2 = obtainnodes(ElemDict,Elset2);

    return intersect(Nset1,Nset2);

end

function obtainnodes(ElemDict,Elset)
#=
The function return the nodes belong to the Elset.
=#
    NodesID = [];
    @inbounds @simd for elemid in Elset
        append!(NodesID,ElemDict[elemid]);
    end
    return unique(NodesID)
end

function obtainnodesinelset(Face_all,Elset)
#=
The function obtainnodesinelset obtain the all faces in the Elset.
=#
    facetotal = size(Face_all,1);
    BlSwitch = zeros(Bool,facetotal);
    @time Threads.@threads for kface = 1:facetotal
        if in(Face_all[kface,1],Elset)
            BlSwitch[kface] = true;
        end
    end
    return findall(BlSwitch)
end