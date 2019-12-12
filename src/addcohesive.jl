function addcohesive_2d!(NodeDict::Dict{Any,Any},ElemDict::Dict{Any,Any},ElsetPartsArray::Array{Array{Int64,1},1},ElsetMatrix::Array{Int64,1})
#=
    The function addcohesive_2d return the new cohesive elements in the dicts and the mapping dictionary between old and new nodes.
=#
    npart = length(ElsetPartsArray);    # obtain total number of element sets
    Face_all,Face_all_Normal = AllFaceGet(NodeDict,ElemDict);
    NodeOld2NewDict = Dict();
    @inbounds @simd for kpart = 1:npart
        addcohesive_2d!(NodeDict,ElemDict,NodeOld2NewDict,ElsetPartsArray[kpart],ElsetMatrix,Face_all=Face_all);
    end
        
    return NodeOld2NewDict
end

function addcohesive_2d!(NodeDict::Dict{Any,Any},ElemDict::Dict{Any,Any},NodeOld2NewDict,ElsetPart::Array{Int64,1},ElsetMatrix::Array{Int64,1}; Face_all = [])
#=
The function addcohesive_2d return the new cohesive elements in the dicts and the mapping dictionary between old and new nodes.
=#  
    
    faceloc = obtainfaceonelset(Face_all[:,1],ElsetPart);
    FaceHere = Face_all[faceloc,:];
    
    # update the elem&node info
    NodebeenModified = collect(NodeOld2NewDict);
    buildcohesive! = buildcohesive!(NodeDict,ElemDict,NodeOld2NewDict,FaceHere);
    NodeModified = setdiff(collect(NodeOld2NewDict),NodebeenModified);

    # replace the modified nodes
    Threads.@threads for elemid in ElsetMatrix
        Nodeshere = ElemDict[elemid];
        @inbounds @simd for knode in 1:length(Nodeshere)
            if Nodeshere[knode] in NodeModified
                ElemDict[elemid][knode] = NodeOld2NewDict[Nodeshere[knode]];
            end
        end
    end

end

function buildcohesive!(NodeDict,ElemDict,NodeOld2NewDict,FaceHere)
#=
  build cohesive element and update the nodedict & elemdict 
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
        NodeOld2NewDict[node2] = node3;

        node4 = nodeidnow;
        nodeidnow += 1;
        NodeDict[node4] = NodeDict[node1];
        NodeOld2NewDict[node1] = node4;

        ElemDict[elemidnow] = [node1;node2;node3;node4];
    end
    
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