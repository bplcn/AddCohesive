function addcohesive_2d!(NodeDict,ElemDict,ElsetPart::Array{Int64,1},ElsetMatrix::Array{Int64,1};Face_all)
#=
The function addcohesive_2d return the new cohesive elements in the dicts and the mapping dictionary between old and new nodes.
=#  
    
    faceloc = obtainfaceonelset(Face_all[:,1],ElsetPart);
    FaceHere = Face_all[faceloc,:];
    buildcohesive! = buildcohesive!(NodeDict,ElemDict,FaceHere);
    
end

function buildcohesive!(NodeDict,ElemDict,FaceHere)
#=
  build cohesive 
=# 
    nodeidnow = maximum(keys(NodeDict))+1;
    elemidnow = maximum(keys(ElemDict))+1;

    nface = size(FaceHere,1)
    for kface = 1:nface
        node1 = FaceHere[kface,3];
        node2 = FaceHere[kface,4];

        node3 = nodeidnow;
        nodeidnow += 1;
        NodeDict[node3] = NodeDict[node2];

        node4 = nodeidnow;
        nodeidnow += 1;
        NodeDict[node4] = NodeDict[node1];

        ElemDict[elemidnow] = [node1;node2;node3;node4];
    end

end

function addcohesive_2d(NodeDict,ElemDict,ElsetPartsArray::Array{Array{Int64,1},1})
#=
    The function addcohesive_2d return the new cohesive elements in the dicts and the mapping dictionary between old and new nodes.
=#
    npart = length(ElsetPartsArray);    # obtain total number of element sets
    Face_all,Face_all_Normal = AllFaceGet(NodeDict,ElemDict);
    
    

    
    

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