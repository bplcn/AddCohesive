using AbaAccess
# obtain the mesh information
include("addcohesive.jl")

InpName = "D:/Abaqus_Temp/temp_near_26nov.inp";
NodeDict,ElemDict,NsetDict,ElsetDict = MeshObtain(InpName);
# got the matrix nodes and elements ids.
NodesIDAll = collect(keys(NodeDict));
ElemsIDAll = collect(keys(ElemDict));

ElemsIDMax = maximum(ElemsIDAll);

# obtain parts and matrix information
ElsetMatrix = deepcopy(ElemsIDAll);
ElsetPartsArray = [];
for kpart = 1:200
    try
        ElsetPart = ElsetDict["CF$(kpart)"];
        push!(ElsetPartsArray,ElsetPart);
        setdiff!(ElsetMatrix,ElsetPart);
    catch
        break
    end
end

using PBCHandler2D
Face_all,Face_all_Normal = AllFaceGet(NodeDict,ElemDict);

addcohesive_2d!(NodeDict,ElemDict,ElsetPartsArray,ElsetMatrix);

InpName = "D:/Abaqus_temp/addcohesive_2d_test.inp"
fID = open(InpName,"w");
# 
NodeIDArray = sort(collect(keys(NodeDict)));
println(fID,"*Node");
for knodeID in NodeIDArray
    println(fID,"$(knodeID),$(NodeDict[knodeID][1]),$(NodeDict[knodeID][2])");
end

ElemIDArray = sort(collect(keys(ElemDict)));
SolNodesIDArray = Array{Int64,1}();
CohNodesIDArray = Array{Int64,1}();
for kelemID in ElemIDArray
    NodeIDHere = ElemDict[kelemID];
    if length(NodeIDHere)==4
        if kelemID>ElemsIDMax 
            push!(CohNodesIDArray,kelemID);
        else
            push!(SolNodesIDArray,kelemID);
        end
    end
    
end
println(fID,"*ELement, type=CPE4, elset=MATRIX")
for kelemID in SolNodesIDArray
    NodeIDHere = ElemDict[kelemID];
    println(fID,"$(kelemID),$(NodeIDHere[1]),$(NodeIDHere[2]),$(NodeIDHere[3]),$(NodeIDHere[4])");
end

println(fID,"*ELement, type=COH2D4, elset=CoH")
for kelemID in CohNodesIDArray
    NodeIDHere = ElemDict[kelemID];
    println(fID,"$(kelemID),$(NodeIDHere[1]),$(NodeIDHere[2]),$(NodeIDHere[3]),$(NodeIDHere[4])");
end

NodeSetNameArray = collect(keys(NsetDict));
for NodeSetName in NodeSetNameArray
    println(fID,"*Nset,nset="*NodeSetName);
    NSetHere = NsetDict[NodeSetName];
    nline = Int64.((length(NSetHere)+10-mod(length(NSetHere),10))/10);
    for kline = 1:(nline-1)
        strtemp = "";
        for kmem = 1:10
            strtemp = strtemp*"$(NSetHere[(kline-1)*10+kmem]),";
        end
        println(fID,strtemp);
    end

    strtemp = "";
    for kmem = 1:mod(length(NSetHere),10)
        strtemp = strtemp*"$(NSetHere[(nline-1)*10+kmem]),";
    end
    if ~isempty(strtemp)
        println(fID,strtemp);
    end
end

ElemSetNameArray = collect(keys(ElsetDict));

for ElemSetName in ElemSetNameArray
    println(fID,"*Elset,elset="*ElemSetName);
    ElSetHere = ElsetDict[ElemSetName];
    nline = Int64.((length(ElSetHere)+10-mod(length(ElSetHere),10))/10);
    for kline = 1:(nline-1)
        strtemp = "";
        for kmem = 1:10
            strtemp = strtemp*"$(ElSetHere[(kline-1)*10+kmem]),";
        end
        println(fID,strtemp);
    end

    strtemp = "";
    for kmem = 1:mod(length(ElSetHere),10)
        strtemp = strtemp*"$(ElSetHere[(nline-1)*10+kmem]),";
    end
    if ~isempty(strtemp)
        println(fID,strtemp);
    end
end
close(fID)
