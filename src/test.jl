using AbaAccess
# obtain the mesh information
InpName = "D:/Abaqus_Temp/temp_near_26nov.inp";
NodeDict,ElemDict,NsetDict,ElsetDict = MeshObtain(InpName);
# got the matrix nodes and elements ids.
NodesIDAll = collect(keys(NodeDict));
ElemsIDAll = collect(keys(ElemDict));
MatrixNodesID = deepcopy(NodesIDAll);
MatrixElemsID = deepcopy(ElemsIDAll);
PartElemsID = [];
for kpart = 1:200
    try
        PartsElemsID = ElsetDict["CF$(kpart)"];
        setdiff!(MatrixElemsID,PartsElemsID);
        append!(PartElemsID,PartsElemsID);
    catch
        break
    end
end

using PBCHandler2D
Face_all,Face_all_Normal = AllFaceGet(NodeDict,ElemDict);
obtainnodesinelset(Face_all[:,1],ElsetDict["CF1"])

faceloc = findall(MatrixSwitch)



# nodesshared = findsharednodes(ElemDict,MatrixElemsID,ElsetDict["CF$(1)"])


