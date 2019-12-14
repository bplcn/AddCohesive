function addcohesive_2d_all!(NodeDict::Dict{Any,Any},ElemDict::Dict{Any,Any},ElsetDict,ElsetPartsArray,ElsetMatrix)
#=
    The function addcohesive_2d return the new cohesive elements in the dicts and the mapping dictionary between old and new nodes.
=#
    npart = length(ElsetPartsArray);    # obtain total number of element sets
    NodesinMatrix = obtainnodes(ElemDict,ElsetMatrix);
    
    Face_all,Face_all_Normal = AllFaceGet(NodeDict,ElemDict);
    NodeOld2NewDict = Dict();
    @inbounds @simd for kpart = 1:npart
        ElemIDOld = collect(keys(ElemDict));
        addcohesive_2d!(NodeDict,ElemDict,NodeOld2NewDict,ElsetPartsArray[kpart],ElsetMatrix,Face_all=Face_all,NodesinMatrix=NodesinMatrix);
        ElemIDNew = collect(keys(ElemDict));
        ElsetDict["Inf_$(kpart)"] = setdiff(ElemIDNew,ElemIDOld);
    end
    
    return NodeOld2NewDict
end

function addcohesive_2d!(NodeDict::Dict{Any,Any},ElemDict::Dict{Any,Any},NodeOld2NewDict,ElsetPart,ElsetMatrix;Face_all,NodesinMatrix)
#=
The function addcohesive_2d return the new cohesive elements in the dicts and the mapping dictionary between old and new nodes.
=#  
    if length(NodesinMatrix) == 0
        NodesinMatrix = obtainnodes(ElemDict,ElsetMatrix);
    end
    
    faceloc = obtainfaceonelset(Face_all[:,1],ElsetPart);
    FaceinPart = Face_all[faceloc,:]
    faceloc = obtainfacehavenodes(FaceinPart,ElemDict,NodesinSet=NodesinMatrix);
    FaceHere = FaceinPart[faceloc,:];
    
    # update the elem&node info
    NodebeenModified = collect(keys(NodeOld2NewDict));
    buildcohesive!(NodeDict,ElemDict,NodeOld2NewDict,FaceHere);
    NodeModified = setdiff(collect(keys(NodeOld2NewDict)),NodebeenModified);

    # replace the modified nodes
    Threads.@threads for elemid in ElsetMatrix
        Nodeshere = copy(ElemDict[elemid]);
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
        if node2 in keys(NodeOld2NewDict)
            node3 = NodeOld2NewDict[node2];
        else
            node3 = nodeidnow;
            nodeidnow += 1;
            NodeDict[node3] = NodeDict[node2];
            NodeOld2NewDict[node2] = node3;
        end

        if node1 in keys(NodeOld2NewDict)
            node4 = NodeOld2NewDict[node1];
        else
            node4 = nodeidnow;
            nodeidnow += 1;
            NodeDict[node4] = NodeDict[node1];
            NodeOld2NewDict[node1] = node4;
        end
        
        ElemDict[elemidnow] = [node1;node2;node3;node4];
        elemidnow +=1
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
        if in(Facehere[kface,3],NodesinSet) && in(Facehere[kface,4],NodesinSet)
            BlSwitch[kface] = true;
        end
    end
    return findall(BlSwitch)

end