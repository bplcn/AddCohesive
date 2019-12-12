using AbaAccess
# obtain the mesh information
include("addcohesive.jl")

InpName = "D:/Abaqus_Temp/temp_near_26nov.inp";
NodeDict,ElemDict,NsetDict,ElsetDict = MeshObtain(InpName);
# got the matrix nodes and elements ids.
NodesIDAll = collect(keys(NodeDict));
ElemsIDAll = collect(keys(ElemDict));

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

# faceloc = obtainnodesinelset(Face_all[:,1],ElsetDict["CF1"]);
# Face_attached_here = Face_all[faceloc,:];

# faceloc = findall(MatrixSwitch)



# nodesshared = findsharednodes(ElemDict,MatrixElemsID,ElsetDict["CF$(1)"])


