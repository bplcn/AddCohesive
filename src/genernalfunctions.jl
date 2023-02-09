
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

function obtainfaceonelset(Face_all,Elset)
#=
    The function obtain the all faces location in Face_all attached to the given Elset.
=#
    facetotal = size(Face_all,1);
    BlSwitch = zeros(Bool,facetotal);
    Threads.@threads for kface = 1:facetotal
        if in(Face_all[kface,1],Elset)
            BlSwitch[kface] = true;
        end
    end
    return findall(BlSwitch)

end

function obtainfacehavenodes(Facehere,ElemDict;Elset)
#=
    The function obtainfacehavenodes return all the face which have nodes in Elset
=#
    NodesinSet = obtainnodes(ElemDict,Elset);
    return obtainfacehavenodes(Facehere,ElemDict,NodesinSet)

end

function obtainfacehavenodes(Facehere,ElemDict;NodesinSet)
#=
    The function obtainfacehavenodes return all the face which have nodes in Elset
=#
        # NodesinSet = obtainnodes(ElemDict,Elset);

    nface = size(Facehere,1);
    BlSwitch = zeros(Bool,nface);
    Threads.@threads for kface in 1:nface
        # if in(Facehere[kface,3],NodesinSet) && in(Facehere[kface,4],NodesinSet)
        #     BlSwitch[kface] = true;
        # end
        BlSwitch[kface] = true
        for knode in 3:length(Facehere[kface,:])
            BlSwitch[kface] = BlSwitch[kface] && in(Facehere[kface,knode],NodesinSet)
        end
    end
    return findall(BlSwitch)

end