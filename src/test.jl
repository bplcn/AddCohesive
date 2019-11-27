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


