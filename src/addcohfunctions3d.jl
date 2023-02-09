function AllFaceGet3d(NodeDict::Dict,ElemDict::Dict)
    `This function get all faces and the normal vector of each elemental triangle surface.`
    # 
    ElemIDArray = collect(keys(ElemDict))
    ELemTotal = length(ElemIDArray)
# 
    Face_1_all = zeros(Int64,ELemTotal,5)
    Face_2_all = zeros(Int64,ELemTotal,5)
    Face_3_all = zeros(Int64,ELemTotal,5)
    Face_4_all = zeros(Int64,ELemTotal,5)

    # @inbounds @simd for kelem = 1:ELemTotal
    Threads.@threads for kelem = 1:ELemTotal
        #
        elemID = ElemIDArray[kelem]

        node1 = ElemDict[elemID][1]
        node2 = ElemDict[elemID][2]
        node3 = ElemDict[elemID][3]
        node4 = ElemDict[elemID][4]

        Face_1_all[kelem,:] = [elemID 1 node1 node3 node2]
        Face_2_all[kelem,:] = [elemID 2 node1 node2 node4]
        Face_3_all[kelem,:] = [elemID 3 node2 node3 node4]
        Face_4_all[kelem,:] = [elemID 4 node1 node4 node3]

        #
    end
    Face_all = [Face_1_all;Face_2_all;Face_3_all;Face_4_all]
    
    return Face_all
end


function addcohesive_3d_all!(NodeDict,ElemDict,ElsetDict,ElsetPartsArray,ElsetMatrix)
#=
    The function addcohesive_2d return the new cohesive elements in the dicts and the mapping dictionary between old and new nodes.
=#  
    npart = length(ElsetPartsArray);    # obtain total number of element sets
    NodesinMatrix = obtainnodes(ElemDict,ElsetMatrix);
    
    Face_all = AllFaceGet3d(NodeDict,ElemDict);
    NodeOld2NewDict = Dict();
    @inbounds @simd for kpart = 1:npart
        ElemIDOld = collect(keys(ElemDict));
        addcohesive_3d!(NodeDict,ElemDict,NodeOld2NewDict,ElsetPartsArray[kpart],ElsetMatrix,Face_all=Face_all,NodesinMatrix=NodesinMatrix);
        ElemIDNew = collect(keys(ElemDict));
        ElsetDict["Inf_$(kpart)"] = setdiff(ElemIDNew,ElemIDOld);
    end
    return NodeOld2NewDict
end

function addcohesive_3d!(NodeDict,ElemDict,NodeOld2NewDict,ElsetPart,ElsetMatrix;Face_all,NodesinMatrix)
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
    buildcohesive3d!(NodeDict,ElemDict,NodeOld2NewDict,FaceHere);
    NodeModified = setdiff(collect(keys(NodeOld2NewDict)),NodebeenModified);

    # replace the modified nodes
    Threads.@threads for elemid in ElsetMatrix
        Nodeshere = copy(ElemDict[elemid]);
        # @inbounds @simd for knode in 1:length(Nodeshere)
        @inbounds @simd for knode in 1:length(Nodeshere)
            if Nodeshere[knode] in NodeModified
                ElemDict[elemid][knode] = NodeOld2NewDict[Nodeshere[knode]];
            end
        end
    end
end

function buildcohesive3d!(NodeDict,ElemDict,NodeOld2NewDict,FaceHere)
#=
    build cohesive element and update the nodedict & elemdict 
=# 
    nodeidnow = maximum(keys(NodeDict))+1;
    elemidnow = maximum(keys(ElemDict))+1;

    nface = size(FaceHere,1)
    for kface = 1:nface
        node1 = FaceHere[kface,3];
        node2 = FaceHere[kface,4];
        node3 = FaceHere[kface,5];
        if node1 in keys(NodeOld2NewDict)
            node4 = NodeOld2NewDict[node1];
        else
            node4 = nodeidnow;
            nodeidnow += 1;
            NodeDict[node4] = NodeDict[node1];
            NodeOld2NewDict[node1] = node4;
        end
        if node2 in keys(NodeOld2NewDict)
            node5 = NodeOld2NewDict[node2];
        else
            node5 = nodeidnow;
            nodeidnow += 1;
            NodeDict[node5] = NodeDict[node2];
            NodeOld2NewDict[node2] = node5;
        end
        if node3 in keys(NodeOld2NewDict)
            node6 = NodeOld2NewDict[node3];
        else
            node6 = nodeidnow;
            nodeidnow += 1;
            NodeDict[node6] = NodeDict[node3];
            NodeOld2NewDict[node3] = node6;
        end
        ElemDict[elemidnow] = [node1;node2;node3;node4;node5;node6];
        elemidnow += 1
    end
end
