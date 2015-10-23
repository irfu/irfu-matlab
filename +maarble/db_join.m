function res = db_join(band,varStr,scList)

res = [];
for iSc=1:length(scList)
  db = evalin('caller',['dbEMIC_' band '_' scList{iSc}]);
  res = [res; db.(varStr)]; %#ok<AGROW>
end